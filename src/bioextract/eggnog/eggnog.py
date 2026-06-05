from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .constant import (
    ASSET_SPECS,
    EggnogInputIdKind,
    MEDIA_TYPE_SQLITE,
    MEDIA_TYPE_SQLITE_GZIP,
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    build_mapping_frame,
    extract_unmapped_input_ids_frame,
    read_cog_fun_frame,
    select_mapping_frame,
    validate_kind_input_id,
)

__all__ = [
    "EggnogDb",
    "EggnogResourceLimits",
    "EggnogSelection",
    "EggnogTidyDataset",
]


@dataclass(frozen=True, slots=True)
class EggnogResourceLimits:
    file_eggnog_db_bytes_max: int | None = None
    file_cog_fun_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _EggnogSnapshot:
    file_eggnog_db: Path
    file_cog_fun: Path | None
    dir_tmp: Path | None


EggnogTidyDataset = TidyDataset


@dataclass(slots=True)
class EggnogDb:
    """Path-first access to local eggNOG mapper resource snapshots."""

    snapshot: _EggnogSnapshot
    limits: EggnogResourceLimits
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_cog_fun: pl.DataFrame | None = field(default=None, init=False, repr=False)

    DEFAULT_RESOURCE_LIMITS = EggnogResourceLimits()

    @classmethod
    def from_files(
        cls,
        *,
        file_eggnog_db: os.PathLike[str] | str,
        file_cog_fun: os.PathLike[str] | str | None = None,
        dir_tmp: os.PathLike[str] | str | None = None,
        limits: EggnogResourceLimits | None = None,
    ) -> EggnogDb:
        """Create a dataset handle from explicit eggNOG mapper files."""
        limits_resolved = EggnogResourceLimits() if limits is None else limits
        file_eggnog_db = _validate_file(
            file_eggnog_db,
            size_max=limits_resolved.file_eggnog_db_bytes_max,
            label="eggNOG SQLite database file",
        )
        if file_cog_fun is not None:
            file_cog_fun = _validate_file(
                file_cog_fun,
                size_max=limits_resolved.file_cog_fun_bytes_max,
                label="COG function table file",
            )
        dir_tmp_resolved = Path(dir_tmp) if dir_tmp is not None else None
        if dir_tmp_resolved is not None:
            dir_tmp_resolved.mkdir(parents=True, exist_ok=True)

        return cls(
            snapshot=_EggnogSnapshot(
                file_eggnog_db=file_eggnog_db,
                file_cog_fun=file_cog_fun,
                dir_tmp=dir_tmp_resolved,
            ),
            limits=limits_resolved,
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full eggNOG protein-to-COG mapping table."""
        if self._df_mapping is None:
            self._df_mapping = build_mapping_frame(
                file_eggnog_db=self.snapshot.file_eggnog_db,
                dir_tmp=self.snapshot.dir_tmp,
                df_cog_fun=self._read_cog_fun(),
            )
        return self._df_mapping

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: EggnogInputIdKind,
    ) -> EggnogSelection:
        """Create a single-query eggNOG mapping selection."""
        validate_kind_input_id(kind_input_id)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return EggnogSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            kind_input_id=kind_input_id,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
        *,
        kind_input_id: EggnogInputIdKind,
    ) -> EggnogSelection:
        """Create a grouped eggNOG mapping selection."""
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
        return EggnogSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            kind_input_id=kind_input_id,
        )

    def build_tidy(self) -> EggnogTidyDataset:
        """Build the in-memory eggNOG tidy dataset."""
        return EggnogTidyDataset(
            frames={"mapping": self.extract_mapping()},
            source=self._tidy_sources(),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"eggnog-mapping-{self.snapshot.file_eggnog_db.stem}",
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
    ) -> TidyWriteReport:
        """Write the eggNOG tidy dataset as flat parquet files."""
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
        )

    def _read_cog_fun(self) -> pl.DataFrame:
        if self._df_cog_fun is None:
            self._df_cog_fun = read_cog_fun_frame(self.snapshot.file_cog_fun)
        return self._df_cog_fun

    def _tidy_sources(self) -> tuple[TidySource, ...]:
        media_type_db = (
            MEDIA_TYPE_SQLITE_GZIP
            if self.snapshot.file_eggnog_db.suffix == ".gz"
            else MEDIA_TYPE_SQLITE
        )
        sources = [
            TidySource(path=self.snapshot.file_eggnog_db, media_type=media_type_db)
        ]
        if self.snapshot.file_cog_fun is not None:
            sources.append(
                TidySource(path=self.snapshot.file_cog_fun, media_type=MEDIA_TYPE_TSV)
            )
        return tuple(sources)


@dataclass(slots=True)
class EggnogSelection:
    """Selection handle for single and grouped eggNOG mapping queries."""

    dataset: EggnogDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: EggnogInputIdKind
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
        """Extract selected eggNOG mapping rows."""
        if self._df_mapping is None:
            self._df_mapping = select_mapping_frame(
                file_eggnog_db=self.dataset.snapshot.file_eggnog_db,
                dir_tmp=self.dataset.snapshot.dir_tmp,
                df_input_ids=self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
                df_cog_fun=self.dataset._read_cog_fun(),
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that did not map to eggNOG rows."""
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
