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
    path: Path
    media_type: str


@dataclass(frozen=True, slots=True)
class TidyAsset:
    path: str
    kind: str
    frame_name: str
    is_optional: bool = False


@dataclass(frozen=True, slots=True)
class TidyReportAsset:
    path: str
    kind: str
    is_optional: bool = False


@dataclass(frozen=True, slots=True)
class TidyManifestAsset:
    path: str
    kind: str
    sha256: str | None
    is_optional: bool = False


class TidyManifest(TypedDict):
    build_id: str
    schema_version: str
    generated_at: str
    sources: list[dict[str, str | int]]
    assets: list[dict[str, str | bool | None]]


@dataclass(frozen=True, slots=True)
class TidyWriteReport:
    dir_out: Path
    assets: tuple[TidyReportAsset, ...]
    manifest: TidyManifest | None = None


@dataclass(slots=True)
class TidyDataset:
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
            self.build_manifest(entries_manifest) if should_write_manifest else None
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

    def build_manifest(
        self,
        assets: list[TidyManifestAsset],
    ) -> TidyManifest:
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
