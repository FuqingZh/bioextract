from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any

import polars as pl

from bioextract._shared import validate_file_size
from bioextract._tidy import TidyWriteReport

from .constant import (
    COLS_IDMAPPING_SELECTED,
    MEDIA_TYPE_PARQUET,
    MEDIA_TYPE_PARQUET_DATASET,
    MEDIA_TYPE_TSV,
    MEDIA_TYPE_TSV_GZIP,
    SCHEMA_VERSION,
)
from .util import (
    filter_taxids,
    has_hive_parquet_candidates,
    normalize_taxids,
    scan_hive_mapping_dataset,
    scan_parquet_mapping,
    scan_raw_idmapping_selected,
    validate_mapping_schema,
)

__all__ = [
    "UniprotDb",
    "UniprotResourceLimits",
]


class _UniprotMappingKind(StrEnum):
    RAW_TSV = "raw_tsv"
    RAW_TSV_GZIP = "raw_tsv_gzip"
    PARQUET = "parquet"
    HIVE_PARQUET = "hive_parquet"


@dataclass(frozen=True, slots=True)
class UniprotResourceLimits:
    file_idmapping_selected_bytes_max: int | None = None


@dataclass(frozen=True, slots=True)
class _UniprotSnapshot:
    file_idmapping_selected: Path
    kind: _UniprotMappingKind
    taxids: tuple[str, ...] = ()


@dataclass(slots=True)
class UniprotDb:
    """Path-first access to UniProt idmapping selected resources.

    `UniprotDb` reads either the raw UniProt `idmapping_selected.tab(.gz)` file,
    a single normalized parquet file, or a hive-partitioned parquet dataset
    written by :meth:`write_tidy`. Construction is deliberately lightweight:
    paths are validated, but data and schemas are not scanned until extraction,
    validation, or writing is requested.
    """

    snapshot: _UniprotSnapshot
    limits: UniprotResourceLimits = field(default_factory=UniprotResourceLimits)

    DEFAULT_RESOURCE_LIMITS = UniprotResourceLimits()

    @classmethod
    def from_files(
        cls,
        *,
        file_idmapping_selected: os.PathLike[str] | str,
        limits: UniprotResourceLimits | None = None,
    ) -> UniprotDb:
        """Create a dataset handle from raw or tidy UniProt mapping data.

        Args:
            file_idmapping_selected: Path to `idmapping_selected.tab.gz`, a
                normalized parquet file, or a hive parquet dataset directory.
            limits: Dataset-level resource limits. Size limits apply to file
                inputs; hive dataset directories are checked structurally only.

        Returns:
            A dataset handle that can be taxid-scoped and materialized later.

        Raises:
            FileNotFoundError: If the path does not exist.
            ValueError: If the path type is unsupported, a directory contains
                no parquet files, or a configured file-size limit is exceeded.
        """
        path = Path(file_idmapping_selected)
        if not path.exists():
            raise FileNotFoundError(
                f"UniProt idmapping selected path not found: {path}"
            )

        limits_resolved = UniprotResourceLimits() if limits is None else limits
        kind = _infer_mapping_kind(path)
        if path.is_file():
            validate_file_size(
                file_path=path,
                size_max=limits_resolved.file_idmapping_selected_bytes_max,
                label="UniProt idmapping selected file",
            )

        return cls(
            snapshot=_UniprotSnapshot(file_idmapping_selected=path, kind=kind),
            limits=limits_resolved,
        )

    def with_taxids(self, *taxids: str | int) -> UniprotDb:
        """Create a taxid-scoped view of this UniProt mapping resource."""
        return UniprotDb(
            snapshot=_UniprotSnapshot(
                file_idmapping_selected=self.snapshot.file_idmapping_selected,
                kind=self.snapshot.kind,
                taxids=normalize_taxids(taxids),
            ),
            limits=self.limits,
        )

    def validate_schema(self) -> None:
        """Validate that the backing data exposes the normalized mapping schema."""
        validate_mapping_schema(self._scan_mapping())

    def extract_mapping(self) -> pl.DataFrame:
        """Extract normalized UniProt idmapping rows for the current taxid scope."""
        lf_mapping = self._scan_mapping()
        validate_mapping_schema(lf_mapping)
        return (
            filter_taxids(lf_mapping, self.snapshot.taxids)
            .select(COLS_IDMAPPING_SELECTED)
            .collect()
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_allow_all: bool = False,
        should_write_manifest: bool = False,
        should_write_hive: bool = True,
        should_overwrite: bool = False,
    ) -> TidyWriteReport:
        """Write normalized UniProt mapping data as parquet.

        Args:
            dir_out: Output directory.
            should_allow_all: Required when no taxids are selected, because all
                taxa may be very large.
            should_write_manifest: Whether to write `manifest.json`.
            should_write_hive: Whether to write hive partition directories.
                This defaults to `True`. Non-hive writing is only supported
                when exactly one taxid is selected.
            should_overwrite: Whether to remove an existing non-empty output
                directory before writing.

        Returns:
            A write report with asset paths and optional manifest content.
        """
        if not self.snapshot.taxids and not should_allow_all:
            raise ValueError(
                "Writing all UniProt taxids requires should_allow_all=True"
            )
        if not should_write_hive and len(self.snapshot.taxids) != 1:
            raise ValueError(
                "should_write_hive=False only supports exactly one selected TaxId"
            )

        dir_out = Path(dir_out)
        if dir_out.exists() and any(dir_out.iterdir()):
            if not should_overwrite:
                raise FileExistsError(
                    f"UniProt tidy output directory is not empty: {dir_out}"
                )
            shutil.rmtree(dir_out)
        dir_out.mkdir(parents=True, exist_ok=True)

        lf_mapping = filter_taxids(self._scan_mapping(), self.snapshot.taxids)
        validate_mapping_schema(lf_mapping)
        if should_write_hive:
            lf_mapping.select(COLS_IDMAPPING_SELECTED).sink_parquet(
                pl.PartitionBy(dir_out, key="TaxId", include_key=True),
                mkdir=True,
            )
            assets = _collect_hive_assets(dir_out)
        else:
            df_mapping = lf_mapping.select(COLS_IDMAPPING_SELECTED).collect()
            file_out = dir_out / "mapping.parquet"
            df_mapping.write_parquet(file_out)
            assets = (
                {
                    "path": "mapping.parquet",
                    "kind": "canonical",
                    "row_count": df_mapping.height,
                    "is_optional": False,
                },
            )

        manifest = self._build_manifest(assets) if should_write_manifest else None
        if manifest is not None:
            (dir_out / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return TidyWriteReport(dir_out=dir_out, assets=assets, manifest=manifest)

    def _scan_mapping(self) -> pl.LazyFrame:
        match self.snapshot.kind:
            case _UniprotMappingKind.RAW_TSV | _UniprotMappingKind.RAW_TSV_GZIP:
                return scan_raw_idmapping_selected(
                    self.snapshot.file_idmapping_selected
                )
            case _UniprotMappingKind.PARQUET:
                return scan_parquet_mapping(self.snapshot.file_idmapping_selected)
            case _UniprotMappingKind.HIVE_PARQUET:
                return scan_hive_mapping_dataset(self.snapshot.file_idmapping_selected)

    def _build_manifest(
        self,
        assets: tuple[dict[str, Any], ...],
    ) -> dict[str, Any]:
        timestamp = datetime.now(UTC)
        source: dict[str, Any] = {
            "path": self.snapshot.file_idmapping_selected.as_posix(),
            "media_type": _media_type_for_kind(self.snapshot.kind),
        }
        if self.snapshot.file_idmapping_selected.is_file():
            source["bytes"] = self.snapshot.file_idmapping_selected.stat().st_size
        return {
            "build_id": f"uniprot-idmapping-selected-{timestamp.strftime('%Y%m%dT%H%M%SZ')}",
            "schema_version": SCHEMA_VERSION,
            "generated_at": timestamp.isoformat().replace("+00:00", "Z"),
            "taxids": list(self.snapshot.taxids),
            "sources": [source],
            "assets": assets,
        }


def _infer_mapping_kind(path: Path) -> _UniprotMappingKind:
    if path.is_dir():
        if not has_hive_parquet_candidates(path):
            raise ValueError(
                f"UniProt hive parquet dataset contains no parquet files: {path}"
            )
        return _UniprotMappingKind.HIVE_PARQUET
    name = path.name
    if name.endswith(".tab.gz"):
        return _UniprotMappingKind.RAW_TSV_GZIP
    if name.endswith(".tab"):
        return _UniprotMappingKind.RAW_TSV
    if path.suffix == ".parquet":
        return _UniprotMappingKind.PARQUET
    raise ValueError(f"Unsupported UniProt idmapping selected input type: {path}")


def _collect_hive_assets(dir_out: Path) -> tuple[dict[str, Any], ...]:
    assets: list[dict[str, Any]] = []
    for file_parquet in sorted(dir_out.rglob("*.parquet")):
        assets.append(
            {
                "path": file_parquet.relative_to(dir_out).as_posix(),
                "kind": "canonical",
                "row_count": None,
                "is_optional": False,
            }
        )
    return tuple(assets)


def _media_type_for_kind(kind: _UniprotMappingKind) -> str:
    match kind:
        case _UniprotMappingKind.RAW_TSV_GZIP:
            return MEDIA_TYPE_TSV_GZIP
        case _UniprotMappingKind.RAW_TSV:
            return MEDIA_TYPE_TSV
        case _UniprotMappingKind.PARQUET:
            return MEDIA_TYPE_PARQUET
        case _UniprotMappingKind.HIVE_PARQUET:
            return MEDIA_TYPE_PARQUET_DATASET
