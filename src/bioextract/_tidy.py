from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import asdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TypedDict

import polars as pl


@dataclass(frozen=True, slots=True)
class TidySource:
    """Describe one source file recorded in a tidy manifest.

    Attributes:
        path: Source file whose path and byte size are recorded.
        media_type: Stable media-type label for the source format.
        sha256: Optional precomputed source digest. When omitted, the manifest
            leaves out the source `sha256` key instead of recording `null`.
    """

    path: Path
    media_type: str
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class TidyAsset:
    """Map one lazy frame to its relative output path and artifact kind.

    `is_optional` is delivery metadata propagated to reports and manifests; it
    does not permit `write()` to skip a missing frame.
    """

    path: str
    kind: str
    frame_name: str
    is_optional: bool = False


@dataclass(frozen=True, slots=True)
class TidyReportAsset:
    """Describe a relative output asset returned after a tidy write."""

    path: str
    kind: str
    is_optional: bool = False


@dataclass(frozen=True, slots=True)
class TidyManifestAsset:
    """Describe an asset in the serializable tidy manifest.

    A missing digest is serialized as a JSON `null` value.
    """

    path: str
    kind: str
    sha256: str | None
    is_optional: bool = False


class TidyManifest(TypedDict):
    """JSON-compatible manifest emitted for one tidy dataset build."""

    build_id: str
    schema_version: str
    generated_at: str
    sources: list[dict[str, str | int]]
    assets: list[dict[str, str | bool | None]]


@dataclass(frozen=True, slots=True)
class TidyWriteReport:
    """Report persisted assets and the optional manifest returned by a write.

    Asset paths are relative to `dir_out`. `manifest` is `None` unless the
    caller requested manifest generation.
    """

    dir_out: Path
    assets: tuple[TidyReportAsset, ...]
    manifest: TidyManifest | None = None


@dataclass(slots=True)
class TidyDataset:
    """Bind lazy resource frames to a versioned flat-file output contract.

    Each asset maps a relative output path to one entry in `frames`. The asset
    tuple defines write and report order; frames not referenced by an asset are
    not persisted. Sources provide provenance for an optional manifest.

    Attributes:
        frames: Lazy frames keyed by the names referenced from `assets`.
        source: One source or an ordered tuple of sources for provenance.
        schema_version: Version of the resource-specific output schema.
        build_id_prefix: Stable prefix used to construct manifest build IDs.
        assets: Ordered output specifications for frames to persist.

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

    def write(
        self,
        dir_out: Path | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Persist configured frames as parquet assets.

        Args:
            dir_out: Output directory. Relative asset paths are resolved below
                this directory.
            should_write_manifest: Whether to also write `manifest.json` and
                return its content.
            should_hash_assets: Whether to calculate SHA256 digests for
                manifest asset entries.

        Returns:
            A report whose assets follow the configured asset order.

        Raises:
            KeyError: If an asset references a frame missing from `frames`.

        Examples:
            Write a one-frame dataset to a temporary directory:

            >>> from tempfile import TemporaryDirectory
            >>> dataset = TidyDataset(
            ...     frames={"term": pl.DataFrame({"id": ["T1"]}).lazy()},
            ...     source=TidySource(Path("data/source.tsv"), "text/tab-separated-values"),
            ...     schema_version="example-v1",
            ...     build_id_prefix="example",
            ...     assets=(TidyAsset("term.parquet", "canonical", "term"),),
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     report = dataset.write(dir_out)
            ...     [asset.path for asset in report.assets]
            ['term.parquet']
        """
        dir_out = Path(dir_out)
        dir_out.mkdir(parents=True, exist_ok=True)

        entries_report: list[TidyReportAsset] = []
        entries_manifest: list[TidyManifestAsset] = []
        for asset in self.assets:
            frame = self.frames[asset.frame_name]
            file_out = dir_out / asset.path
            file_out.parent.mkdir(parents=True, exist_ok=True)
            frame.sink_parquet(file_out)
            entry_report = TidyReportAsset(
                path=asset.path,
                kind=asset.kind,
                is_optional=asset.is_optional,
            )
            entries_report.append(entry_report)
            entries_manifest.append(
                TidyManifestAsset(
                    path=asset.path,
                    kind=asset.kind,
                    is_optional=asset.is_optional,
                    sha256=calculate_file_sha256(file_out)
                    if should_hash_assets
                    else None,
                )
            )

        manifest = (
            self._build_manifest(entries_manifest) if should_write_manifest else None
        )
        if manifest is not None:
            (dir_out / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return TidyWriteReport(
            dir_out=dir_out,
            assets=tuple(entries_report),
            manifest=manifest,
        )

    def _build_manifest(
        self,
        assets: list[TidyManifestAsset],
    ) -> TidyManifest:
        """Build manifest content from source metadata and prepared assets.

        Args:
            assets: Ordered asset records, including any digests calculated at
                the write boundary.

        Returns:
            JSON-compatible manifest content with a UTC build timestamp.

        Notes:
            Manifest construction stays behind `write()` so callers cannot
            bypass asset persistence or request hashes that were never
            calculated. This method does not read assets to calculate hashes.
        """
        timestamp = datetime.now(UTC)
        return {
            "build_id": f"{self.build_id_prefix}-{timestamp.strftime('%Y%m%dT%H%M%SZ')}",
            "schema_version": self.schema_version,
            "generated_at": timestamp.isoformat().replace("+00:00", "Z"),
            "sources": [
                {
                    "path": source.path.as_posix(),
                    "bytes": source.path.stat().st_size,
                    "media_type": source.media_type,
                    **({"sha256": source.sha256} if source.sha256 is not None else {}),
                }
                for source in self._sources
            ],
            "assets": [asdict(asset) for asset in assets],
        }

    @property
    def _sources(self) -> tuple[TidySource, ...]:
        if isinstance(self.source, TidySource):
            return (self.source,)
        return self.source


def calculate_file_sha256(file_path: Path) -> str:
    with file_path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()
