from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .brite.constant import (
    ASSET_SPECS as BRITE_ASSET_SPECS,
    MEDIA_TYPE_JSON,
    SCHEMA_VERSION as BRITE_SCHEMA_VERSION,
)
from .brite.tidy import build_tidy_frames as build_brite_tidy_frames
from .mapping.constant import (
    ASSET_SPECS as MAPPING_ASSET_SPECS,
    KeggInputIdKind,
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION as MAPPING_SCHEMA_VERSION,
)
from .mapping.util import (
    build_mapping_frame,
    extract_mapping_frame,
    extract_unmapped_input_ids_frame,
    read_conv_ncbi_geneid_frame,
    read_conv_uniprot_frame,
    read_gene_ko_frame,
    read_gene_list_frame,
    read_gene_pathway_frame,
    validate_kind_input_id,
)

__all__ = [
    "KeggDb",
    "KeggResourceLimits",
    "KeggTidyDataset",
]


class _KeggSnapshotKind(StrEnum):
    BRITE_JSON = "brite_json"
    MAPPING_FILES = "mapping_files"


@dataclass(frozen=True, slots=True)
class KeggResourceLimits:
    file_brite_json_bytes_max: int | None = None
    file_conv_uniprot_bytes_max: int | None = None
    file_gene_ko_bytes_max: int | None = None
    file_gene_pathway_bytes_max: int | None = None
    file_gene_list_bytes_max: int | None = None
    file_conv_ncbi_geneid_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _KeggSnapshot:
    kind: _KeggSnapshotKind
    file_brite_json: Path | None = None
    file_conv_uniprot: Path | None = None
    file_gene_ko: Path | None = None
    file_gene_pathway: Path | None = None
    organism_code: str | None = None
    file_gene_list: Path | None = None
    file_conv_ncbi_geneid: Path | None = None


KeggTidyDataset = TidyDataset


@dataclass(slots=True)
class KeggDb:
    """Path-first access to local KEGG resource snapshots."""

    snapshot: _KeggSnapshot
    limits: KeggResourceLimits
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)

    DEFAULT_RESOURCE_LIMITS = KeggResourceLimits()

    @classmethod
    def from_brite_json(
        cls,
        file_brite_json: os.PathLike[str] | str,
        *,
        limits: KeggResourceLimits | None = None,
    ) -> KeggDb:
        """Create a dataset handle from a local KEGG BRITE JSON file."""
        limits_resolved = KeggResourceLimits() if limits is None else limits
        file_brite_json = _validate_file(
            file_brite_json,
            size_max=limits_resolved.file_brite_json_bytes_max,
            label="KEGG BRITE JSON file",
        )
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.BRITE_JSON,
                file_brite_json=file_brite_json,
            ),
            limits=limits_resolved,
        )

    @classmethod
    def from_mapping_files(
        cls,
        *,
        file_conv_uniprot: os.PathLike[str] | str,
        file_gene_ko: os.PathLike[str] | str,
        file_gene_pathway: os.PathLike[str] | str,
        organism_code: str,
        file_gene_list: os.PathLike[str] | str | None = None,
        file_conv_ncbi_geneid: os.PathLike[str] | str | None = None,
        limits: KeggResourceLimits | None = None,
    ) -> KeggDb:
        """Create a dataset handle from explicit KEGG organism mapping files."""
        organism_code = str(organism_code).strip()
        if not organism_code:
            raise ValueError("KEGG organism_code must be non-empty after normalization")

        limits_resolved = KeggResourceLimits() if limits is None else limits
        file_conv_uniprot = _validate_file(
            file_conv_uniprot,
            size_max=limits_resolved.file_conv_uniprot_bytes_max,
            label="KEGG conv_uniprot file",
        )
        file_gene_ko = _validate_file(
            file_gene_ko,
            size_max=limits_resolved.file_gene_ko_bytes_max,
            label="KEGG gene_ko file",
        )
        file_gene_pathway = _validate_file(
            file_gene_pathway,
            size_max=limits_resolved.file_gene_pathway_bytes_max,
            label="KEGG gene_pathway file",
        )
        if file_gene_list is not None:
            file_gene_list = _validate_file(
                file_gene_list,
                size_max=limits_resolved.file_gene_list_bytes_max,
                label="KEGG gene_list file",
            )
        if file_conv_ncbi_geneid is not None:
            file_conv_ncbi_geneid = _validate_file(
                file_conv_ncbi_geneid,
                size_max=limits_resolved.file_conv_ncbi_geneid_bytes_max,
                label="KEGG conv_ncbi_geneid file",
            )

        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.MAPPING_FILES,
                file_conv_uniprot=file_conv_uniprot,
                file_gene_ko=file_gene_ko,
                file_gene_pathway=file_gene_pathway,
                organism_code=organism_code,
                file_gene_list=file_gene_list,
                file_conv_ncbi_geneid=file_conv_ncbi_geneid,
            ),
            limits=limits_resolved,
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full KEGG organism mapping table."""
        self._require_mapping_snapshot("extract KEGG mapping")
        if self._df_mapping is None:
            self._df_mapping = build_mapping_frame(
                organism_code=self.snapshot.organism_code or "",
                df_conv_uniprot=read_conv_uniprot_frame(
                    self._required_path(self.snapshot.file_conv_uniprot)
                ),
                df_conv_ncbi_geneid=read_conv_ncbi_geneid_frame(
                    self.snapshot.file_conv_ncbi_geneid
                ),
                df_gene_ko=read_gene_ko_frame(
                    self._required_path(self.snapshot.file_gene_ko)
                ),
                df_gene_pathway=read_gene_pathway_frame(
                    self._required_path(self.snapshot.file_gene_pathway)
                ),
                df_gene_list=read_gene_list_frame(self.snapshot.file_gene_list),
            )
        return self._df_mapping

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: KeggInputIdKind,
    ) -> KeggSelection:
        """Create a single-query KEGG mapping selection."""
        self._require_mapping_snapshot("select KEGG IDs")
        validate_kind_input_id(kind_input_id)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return KeggSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            kind_input_id=kind_input_id,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
        *,
        kind_input_id: KeggInputIdKind,
    ) -> KeggSelection:
        """Create a grouped KEGG mapping selection."""
        self._require_mapping_snapshot("select grouped KEGG IDs")
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
        return KeggSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            kind_input_id=kind_input_id,
        )

    def build_tidy(self) -> KeggTidyDataset:
        """Build the lazy KEGG tidy dataset for this snapshot kind."""
        if self.snapshot.kind == _KeggSnapshotKind.BRITE_JSON:
            file_brite_json = self._required_path(self.snapshot.file_brite_json)
            frames = {
                frame_name: frame.lazy()
                for frame_name, frame in build_brite_tidy_frames(
                    file_brite_json
                ).items()
            }
            return KeggTidyDataset(
                frames=frames,
                source=TidySource(path=file_brite_json, media_type=MEDIA_TYPE_JSON),
                schema_version=BRITE_SCHEMA_VERSION,
                build_id_prefix=f"kegg-brite-{file_brite_json.stem}",
                assets=tuple(
                    TidyAsset(path=path, kind=kind, frame_name=frame_name)
                    for path, kind, frame_name in BRITE_ASSET_SPECS
                ),
            )

        self._require_mapping_snapshot("build KEGG mapping tidy dataset")
        return KeggTidyDataset(
            frames={"mapping": self.extract_mapping().lazy()},
            source=self._mapping_tidy_sources(),
            schema_version=MAPPING_SCHEMA_VERSION,
            build_id_prefix=f"kegg-mapping-{self.snapshot.organism_code}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in MAPPING_ASSET_SPECS
            ),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the KEGG tidy dataset as flat parquet files."""
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def _mapping_tidy_sources(self) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                path=self._required_path(self.snapshot.file_conv_uniprot),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                path=self._required_path(self.snapshot.file_gene_ko),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                path=self._required_path(self.snapshot.file_gene_pathway),
                media_type=MEDIA_TYPE_TSV,
            ),
        ]
        if self.snapshot.file_gene_list is not None:
            sources.append(
                TidySource(path=self.snapshot.file_gene_list, media_type=MEDIA_TYPE_TSV)
            )
        if self.snapshot.file_conv_ncbi_geneid is not None:
            sources.append(
                TidySource(
                    path=self.snapshot.file_conv_ncbi_geneid,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        return tuple(sources)

    def _require_mapping_snapshot(self, action: str) -> None:
        if self.snapshot.kind != _KeggSnapshotKind.MAPPING_FILES:
            raise ValueError(f"Cannot {action} from a KEGG BRITE JSON snapshot")

    @staticmethod
    def _required_path(path: Path | None) -> Path:
        if path is None:
            raise ValueError("Required KEGG resource path is missing")
        return path


@dataclass(slots=True)
class KeggSelection:
    """Selection handle for single and grouped KEGG mapping queries."""

    dataset: KeggDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: KeggInputIdKind
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
        """Extract selected KEGG mapping rows."""
        if self._df_mapping is None:
            self._df_mapping = extract_mapping_frame(
                self.dataset.extract_mapping(),
                self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that did not map to KEGG rows."""
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
