from __future__ import annotations

import hashlib
import os
from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    write_duckdb_publication,
)


@dataclass(frozen=True, slots=True)
class TidySource:
    """Describe one source file recorded in publication provenance.

    Attributes:
        logical_name: Required role name that identifies this source within the
            publication independently of its path or basename.
        path: Source file whose path and byte size are recorded.
        media_type: Stable media-type label for the source format.
        sha256: Optional precomputed source digest.
    """

    logical_name: str
    path: Path
    media_type: str
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class TidyAsset:
    """Map one lazy frame to its resource-local relation metadata.

    `path` remains an internal relation label during the 1.0 convergence; it
    never determines a caller output path. `kind` becomes the DuckDB relation
    role when a multi-relation publication is built.
    """

    path: str
    kind: str
    frame_name: str
    is_optional: bool = False


@dataclass(slots=True)
class TidyDataset:
    """Bind lazy resource relations to versioned publication contracts.

    Each asset maps a relative output path to one entry in `frames`. The asset
    tuple defines write and report order; frames not referenced by an asset are
    not persisted. Sources provide embedded publication provenance.

    Attributes:
        frames: Lazy frames keyed by the names referenced from `assets`.
        source: One source or an ordered tuple of sources for provenance.
        resource_schema_version: Version of the resource-specific output schema.
        source_schema_profile: Required bioextract-owned, content-validated
            parser profile.
        build_id_prefix: Human-readable build/execution label prefix. It never
            supplies resource, release, or schema identity.
        assets: Ordered output specifications for frames to persist.
        resource_name: Stable database resource name for embedded provenance.
        release_version: Optional official source release identifier.
        release_version_source: Optional `caller` or `official_metadata`
            provenance for `release_version`.
        source_schema_version: Optional upstream-declared input schema label.

    Examples:
        Bind one in-memory frame to a parquet asset:

        >>> dataset = TidyDataset(
        ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
        ...     source=TidySource(
        ...         "source", Path("data/source.tsv"), "text/tab-separated-values"
        ...     ),
        ...     resource_schema_version="example-v1",
        ...     source_schema_profile="example-source-v1",
        ...     build_id_prefix="example",
        ...     assets=(TidyAsset("term.parquet", "canonical", "term"),),
        ... )
        >>> dataset.frames["term"].collect().to_dicts()
        [{'id': 'T1'}]
    """

    frames: Mapping[str, pl.LazyFrame]
    source: TidySource | tuple[TidySource, ...]
    resource_schema_version: str
    source_schema_profile: str
    build_id_prefix: str
    assets: tuple[TidyAsset, ...]
    resource_name: str | None = None
    release_version: str | None = None
    release_version_source: str | None = None
    source_schema_version: str | None = None
    scope: str | None = None
    extra_metadata: Mapping[str, str] | None = None
    before_duckdb_commit: Callable[[], None] | None = None

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        table_names: Mapping[str, str] | None = None,
        if_exists: str = "fail",
        source_columns: Mapping[str, Collection[str]] | None = None,
        include_source_hashes: bool = False,
    ) -> DuckDBWriteResult:
        """Publish all configured relations as one provenance-aware DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> dataset = TidyDataset(
            ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
            ...     source=TidySource(
            ...         "source", Path("data/source.tsv"), "text/plain"
            ...     ),
            ...     resource_schema_version="example-v1",
            ...     source_schema_profile="example-source-v1",
            ...     build_id_prefix="example",
            ...     assets=(TidyAsset("term.parquet", "canonical", "term"),),
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     dataset.write_duckdb(
            ...         Path(dir_out) / "example.duckdb"
            ...     ).tables
            ('term',)
        """
        names = table_names or {}
        declared_source_columns = source_columns or {}
        relations = tuple(
            RelationSpec(
                table_name=names.get(asset.frame_name, asset.frame_name),
                frame=self.frames[asset.frame_name],
                role=asset.kind,
                source_columns=tuple(declared_source_columns.get(asset.frame_name, ())),
            )
            for asset in self.assets
        )
        return write_duckdb_publication(
            relations,
            path,
            resource_name=self.resource_name or self.build_id_prefix,
            resource_schema_version=self.resource_schema_version,
            source_schema_profile=self.source_schema_profile,
            source_schema_version=self.source_schema_version,
            sources=self._source_records_with_hashes(include_source_hashes),
            scope=self.scope,
            release_version=self.release_version,
            release_version_source=self.release_version_source,
            if_exists=if_exists,
            extra_metadata=self.extra_metadata,
            before_commit=self.before_duckdb_commit,
        )

    @property
    def _sources(self) -> tuple[TidySource, ...]:
        if isinstance(self.source, TidySource):
            return (self.source,)
        return self.source

    @property
    def _source_records(self) -> tuple[SourceFileRecord, ...]:
        records: list[SourceFileRecord] = []
        logical_names: set[str] = set()
        for source in self._sources:
            logical_name = source.logical_name.strip()
            if not logical_name:
                raise ValueError("TidySource.logical_name must be non-empty")
            if logical_name in logical_names:
                raise ValueError(
                    "TidySource.logical_name values must be unique after normalization"
                )
            logical_names.add(logical_name)
            records.append(
                SourceFileRecord(
                    logical_name=logical_name,
                    path=source.path,
                    media_type=source.media_type,
                    sha256=source.sha256,
                )
            )
        return tuple(records)

    def _source_records_with_hashes(
        self, include_source_hashes: bool
    ) -> tuple[SourceFileRecord, ...]:
        records = self._source_records
        if not include_source_hashes:
            return records
        return tuple(
            SourceFileRecord(
                logical_name=record.logical_name,
                path=record.path,
                media_type=record.media_type,
                sha256=record.sha256 or calculate_file_sha256(record.path),
            )
            for record in records
        )


def calculate_file_sha256(file_path: Path) -> str:
    with file_path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()
