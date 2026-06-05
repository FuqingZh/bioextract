from __future__ import annotations

import json
import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import (
    TidyAsset,
    TidyDataset,
    TidyManifest,
    TidyReportAsset,
    TidySource,
    TidyWriteReport,
    calculate_file_sha256,
)

from .constant import (
    ASSET_SPECS,
    InterProInputIdKind,
    MEDIA_TYPE_TSV_GZIP,
    MEDIA_TYPE_XML_GZIP,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    extract_unmapped_input_ids_frame,
    read_interpro_xml_frames,
    read_mapping_frame,
    scan_mapping_frame,
    select_mapping_frame,
    validate_kind_input_id,
)

__all__ = [
    "InterProDb",
    "InterProResourceLimits",
    "InterProSelection",
    "InterProTidyDataset",
]


@dataclass(frozen=True, slots=True)
class InterProResourceLimits:
    file_protein2ipr_bytes_max: int | None = None
    file_interpro_xml_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _InterProSnapshot:
    file_protein2ipr: Path
    file_interpro_xml: Path | None = None


InterProTidyDataset = TidyDataset


@dataclass(slots=True)
class InterProDb:
    """Path-first access to local InterPro protein mapping snapshots."""

    snapshot: _InterProSnapshot
    limits: InterProResourceLimits = field(default_factory=InterProResourceLimits)
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _frames_xml: dict[str, pl.DataFrame] | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_RESOURCE_LIMITS = InterProResourceLimits()

    @classmethod
    def from_mapping_files(
        cls,
        *,
        file_protein2ipr: os.PathLike[str] | str,
        file_interpro_xml: os.PathLike[str] | str | None = None,
        limits: InterProResourceLimits | None = None,
    ) -> InterProDb:
        """Create a dataset handle from local InterPro mapping files."""
        limits_resolved = InterProResourceLimits() if limits is None else limits
        file_protein2ipr = _validate_file(
            file_protein2ipr,
            size_max=limits_resolved.file_protein2ipr_bytes_max,
            label="InterPro protein2ipr file",
        )
        if file_interpro_xml is not None:
            file_interpro_xml = _validate_file(
                file_interpro_xml,
                size_max=limits_resolved.file_interpro_xml_bytes_max,
                label="InterPro XML file",
            )
        return cls(
            snapshot=_InterProSnapshot(
                file_protein2ipr=file_protein2ipr,
                file_interpro_xml=file_interpro_xml,
            ),
            limits=limits_resolved,
        )

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: InterProInputIdKind,
    ) -> InterProSelection:
        """Create a single-query selection from UniProt accessions."""
        validate_kind_input_id(kind_input_id)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return InterProSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            kind_input_id=kind_input_id,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
        *,
        kind_input_id: InterProInputIdKind,
    ) -> InterProSelection:
        """Create a grouped selection from multiple UniProt accession sets."""
        validate_kind_input_id(kind_input_id)
        grp_in_frames = create_group_input_frames(
            group_to_ids,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        validate_count_limit(
            count=grp_in_frames.df_groups.height,
            limit_max=self.limits.num_groups_max,
            label="Group count",
        )
        validate_count_limit(
            count=grp_in_frames.df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return InterProSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            kind_input_id=kind_input_id,
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full UniProt-to-InterPro mapping table."""
        if self._df_mapping is None:
            self._df_mapping = read_mapping_frame(
                self.snapshot.file_protein2ipr,
                df_interpro_entry=self._xml_frame("entry"),
                df_interpro_member=self._xml_frame("member"),
            )
        return self._df_mapping

    def build_tidy(self) -> InterProTidyDataset:
        """Build the lazy InterPro tidy dataset."""
        return InterProTidyDataset(
            frames={
                "mapping": scan_mapping_frame(
                    self.snapshot.file_protein2ipr,
                    df_interpro_entry=self._xml_frame("entry"),
                    df_interpro_member=self._xml_frame("member"),
                )
            },
            source=self._tidy_sources(),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"interpro-mapping-{self.snapshot.file_protein2ipr.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the InterPro tidy dataset as a flat parquet file."""
        dir_out = Path(dir_out)
        dir_out.mkdir(parents=True, exist_ok=True)
        file_out = dir_out / "mapping.parquet"
        scan_mapping_frame(
            self.snapshot.file_protein2ipr,
            df_interpro_entry=self._xml_frame("entry"),
            df_interpro_member=self._xml_frame("member"),
        ).sink_parquet(file_out)

        asset: TidyReportAsset = {
            "path": "mapping.parquet",
            "kind": "canonical",
            "is_optional": False,
        }
        manifest = (
            self._build_manifest(
                {
                    **asset,
                    "sha256": calculate_file_sha256(file_out)
                    if should_hash_assets
                    else None,
                }
            )
            if should_write_manifest
            else None
        )
        if manifest is not None:
            (dir_out / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return TidyWriteReport(dir_out=dir_out, assets=(asset,), manifest=manifest)

    def _build_manifest(
        self,
        asset: dict[str, str | bool | None],
    ) -> TidyManifest:
        timestamp = datetime.now(UTC)
        return {
            "build_id": f"interpro-mapping-{timestamp.strftime('%Y%m%dT%H%M%SZ')}",
            "schema_version": SCHEMA_VERSION,
            "generated_at": timestamp.isoformat().replace("+00:00", "Z"),
            "sources": [
                {
                    "path": source.path.as_posix(),
                    "bytes": source.path.stat().st_size,
                    "media_type": source.media_type,
                }
                for source in self._tidy_sources()
            ],
            "assets": [asset],
        }

    def _xml_frame(self, frame_name: str) -> pl.DataFrame:
        if self._frames_xml is None:
            self._frames_xml = read_interpro_xml_frames(self.snapshot.file_interpro_xml)
        return self._frames_xml[frame_name]

    def _tidy_sources(self) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                path=self.snapshot.file_protein2ipr, media_type=MEDIA_TYPE_TSV_GZIP
            )
        ]
        if self.snapshot.file_interpro_xml is not None:
            sources.append(
                TidySource(
                    path=self.snapshot.file_interpro_xml, media_type=MEDIA_TYPE_XML_GZIP
                )
            )
        return tuple(sources)


@dataclass(slots=True)
class InterProSelection:
    """Selection handle for single and grouped InterPro mapping queries."""

    dataset: InterProDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: InterProInputIdKind
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs."""
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract selected UniProt-to-InterPro mapping rows."""
        if self._df_mapping is None:
            self._df_mapping = select_mapping_frame(
                self.dataset.snapshot.file_protein2ipr,
                self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
                df_interpro_entry=self.dataset._xml_frame("entry"),
                df_interpro_member=self.dataset._xml_frame("member"),
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that did not map to InterPro rows."""
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmapped_input_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                cols_group_id=self._col_group_id,
            )
        return self._df_unmapped


def _validate_file(
    file_path: os.PathLike[str] | str,
    *,
    size_max: int | None,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    validate_file_size(file_path=file_path, size_max=size_max, label=label)
    return file_path
