from __future__ import annotations

import hashlib
import os
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    ParquetWriteResult,
    RelationSpec,
    SourceFileRecord,
    write_duckdb_publication,
    write_parquet_publication,
)


@dataclass(frozen=True, slots=True)
class TidySource:
    """Describe one source file recorded in publication provenance.

    Attributes:
        path: Source file whose path and byte size are recorded.
        media_type: Stable media-type label for the source format.
        sha256: Optional precomputed source digest.
    """

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
        schema_version: Version of the resource-specific output schema.
        build_id_prefix: Stable resource identity prefix.
        assets: Ordered output specifications for frames to persist.
        resource_name: Stable database resource name for embedded provenance.
        release_version: Optional official source release identifier.

    Examples:
        Bind one in-memory frame to a parquet asset:

        >>> dataset = TidyDataset(
        ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
        ...     source=TidySource(Path("data/source.tsv"), "text/tab-separated-values"),
        ...     schema_version="example-v1",
        ...     build_id_prefix="example",
        ...     assets=(TidyAsset("term.parquet", "canonical", "term"),),
        ... )
        >>> dataset.frames["term"].collect().to_dicts()
        [{'id': 'T1'}]
    """

    frames: Mapping[str, pl.LazyFrame]
    source: TidySource | tuple[TidySource, ...]
    schema_version: str
    build_id_prefix: str
    assets: tuple[TidyAsset, ...]
    resource_name: str | None = None
    release_version: str | None = None

    def write_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
        preserve_source_headers: bool = False,
    ) -> ParquetWriteResult:
        """Publish a single-relation dataset as one provenance-aware Parquet.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> dataset = TidyDataset(
            ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
            ...     source=TidySource(Path("data/source.tsv"), "text/plain"),
            ...     schema_version="example-v1",
            ...     build_id_prefix="example",
            ...     assets=(TidyAsset("term.parquet", "canonical", "term"),),
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     dataset.write_parquet(
            ...         Path(dir_out) / "example.parquet"
            ...     ).schema_version
            'example-v1'
        """
        if len(self.assets) != 1:
            raise ValueError(
                "write_parquet() requires exactly one published relation; "
                "use write_duckdb() for a multi-relation dataset"
            )
        asset = self.assets[0]
        return write_parquet_publication(
            self.frames[asset.frame_name],
            path,
            resource_name=self.resource_name or self.build_id_prefix,
            schema_version=self.schema_version,
            sources=self._source_records,
            release_version=self.release_version,
            if_exists=if_exists,
            normalize_columns=not preserve_source_headers,
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        table_names: Mapping[str, str] | None = None,
        if_exists: str = "fail",
        preserve_source_headers: Collection[str] = (),
    ) -> DuckDBWriteResult:
        """Publish all configured relations as one provenance-aware DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> dataset = TidyDataset(
            ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
            ...     source=TidySource(Path("data/source.tsv"), "text/plain"),
            ...     schema_version="example-v1",
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
        relations = tuple(
            RelationSpec(
                table_name=names.get(asset.frame_name, asset.frame_name),
                frame=self.frames[asset.frame_name],
                role=asset.kind,
                preserve_source_headers=asset.frame_name in preserve_source_headers,
            )
            for asset in self.assets
        )
        return write_duckdb_publication(
            relations,
            path,
            resource_name=self.resource_name or self.build_id_prefix,
            schema_version=self.schema_version,
            sources=self._source_records,
            release_version=self.release_version,
            if_exists=if_exists,
        )

    @property
    def _sources(self) -> tuple[TidySource, ...]:
        if isinstance(self.source, TidySource):
            return (self.source,)
        return self.source

    @property
    def _source_records(self) -> tuple[SourceFileRecord, ...]:
        return tuple(
            SourceFileRecord(
                logical_name=source.path.name,
                path=source.path,
                media_type=source.media_type,
                sha256=source.sha256,
            )
            for source in self._sources
        )


def calculate_file_sha256(file_path: Path) -> str:
    with file_path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()
