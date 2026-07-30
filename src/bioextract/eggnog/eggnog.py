from __future__ import annotations

import os
import tempfile
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
    MEDIA_TYPE_SQLITE,
    MEDIA_TYPE_SQLITE_GZIP,
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
    EggnogInputIdKind,
)
from .util import (
    build_mapping_frame,
    extract_unmapped_input_ids_frame,
    read_cog_fun_frame,
    scan_mapping_tsv,
    select_mapping_frame,
    validate_kind_input_id,
    write_mapping_tsv,
)

__all__ = [
    "EggnogDb",
    "EggnogResourceLimits",
    "EggnogSelection",
    "EggnogTidyDataset",
]


@dataclass(frozen=True, slots=True)
class EggnogResourceLimits:
    """Optional guards for eggNOG source files and query cardinality.

    Attributes:
        file_eggnog_db_bytes_max: Maximum on-disk size of the SQLite or
            gzip-wrapped SQLite source, in bytes. `None` disables the limit.
        file_cog_fun_bytes_max: Maximum on-disk size of the optional COG
            function table, in bytes. `None` disables the limit.
        num_input_ids_max: Maximum number of distinct normalized IDs in one
            selection. `None` disables the limit.
        num_groups_max: Maximum number of groups in one grouped selection.
            `None` disables the limit.

    Examples:
        Reject an oversized eggNOG snapshot before opening SQLite:

        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_db = Path(dir_tmp) / "eggnog.db"
        ...     _ = file_db.write_bytes(b"not sqlite")
        ...     limits = EggnogResourceLimits(file_eggnog_db_bytes_max=1)
        ...     try:
        ...         EggnogDb.from_files(file_eggnog_db=file_db, limits=limits)
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

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
    """Access one local eggNOG mapper resource snapshot.

    Construction validates paths and configured size limits without expanding
    the SQLite mapping. The optional COG function table enriches `CogClass` and
    `CogName`; those columns remain null when the table is omitted. Materialized
    mapping and lookup frames are cached on the handle.

    Examples:
        Read one enriched protein-to-COG annotation:

        >>> db = EggnogDb.from_files(
        ...     file_eggnog_db="fixtures/eggnog/eggnog.db",
        ...     file_cog_fun="fixtures/eggnog/cog-24.fun.tab",
        ... )
        >>> db.extract_mapping().select(
        ...     "EggnogProteinId", "CogCategory", "CogName"
        ... ).head(1).rows()
        [('9606.ENSP1', 'E', 'Amino acid transport and metabolism')]
    """

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
        """Create a handle from explicit eggNOG mapper files.

        Args:
            file_eggnog_db: Path to an eggNOG mapper SQLite database or its
                `.gz` wrapper.
            file_cog_fun: Optional COG function table used to populate class
                and display-name columns.
            dir_tmp: Optional scratch directory for decompressing a wrapped
                SQLite database and staging full tidy exports.
            limits: Optional file-size and selection-cardinality guards.

        Returns:
            A lightweight handle that defers database reads until extraction.

        Raises:
            FileNotFoundError: If a supplied source file does not exist.
            ValueError: If a supplied file exceeds its configured size limit.

        Examples:
            Enrich COG categories with names from the optional lookup:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db",
            ...     file_cog_fun="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.extract_mapping().select(
            ...     "CogCategory", "CogName"
            ... ).head(1).rows()
            [('E', 'Amino acid transport and metabolism')]
        """
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
        """Extract the full protein-to-COG mapping table.

        Returns:
            A cached frame that preserves protein-to-OG and OG-to-COG
            many-to-many relationships. Missing optional COG metadata remains
            null.

        Examples:
            Preserve one row for each protein, OG, and COG-category match:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db",
            ...     file_cog_fun="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.extract_mapping().select(
            ...     "EggnogProteinId", "EggnogOgId", "CogCategory"
            ... ).rows()
            [('9606.ENSP1', 'OG0001', 'E'), ('9606.ENSP1', 'OG0001', 'G'), ('9606.ENSP1', 'OG0002', 'S')]
        """
        if self._df_mapping is None:
            self._df_mapping = build_mapping_frame(
                file_eggnog_db=self.snapshot.file_eggnog_db,
                dir_tmp=self.snapshot.dir_tmp,
                df_cog_fun=self.read_cog_fun(),
            )
        return self._df_mapping

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: EggnogInputIdKind,
    ) -> EggnogSelection:
        """Create a selection for one normalized ID set.

        Args:
            ids: eggNOG protein IDs. Empty values are discarded and duplicate
                normalized IDs collapse to one input row.
            kind_input_id: Input namespace. The supported value is
                `"eggnog_protein"`.

        Returns:
            A selection handle whose outputs include `InputId` and
            `KindInputId` provenance columns.

        Raises:
            ValueError: If the namespace is unsupported or the normalized ID
                count exceeds the configured limit.

        Examples:
            Return annotations only for the selected protein ID:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(
            ...     ["9606.ENSP1"], kind_input_id="eggnog_protein"
            ... )
            >>> selection.extract_mapping().select(
            ...     "InputId", "EggnogOgId", "CogCategory"
            ... ).rows()
            [('9606.ENSP1', 'OG0001', 'E'), ('9606.ENSP1', 'OG0001', 'G'), ('9606.ENSP1', 'OG0002', 'S')]
        """
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
        """Create a selection that preserves caller-defined groups.

        Args:
            group_to_ids: Mapping from unique, non-empty group IDs to eggNOG
                protein IDs.
            kind_input_id: Input namespace. The supported value is
                `"eggnog_protein"`.

        Returns:
            A selection handle whose mapping and unmapped outputs retain
            `GroupId`.

        Raises:
            ValueError: If group IDs or the namespace are invalid, or a group
                or input-count limit is exceeded.

        Examples:
            Preserve a comparison label on every selected annotation row:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["9606.ENSP1"], "down": ["9606.MISSING"]},
            ...     kind_input_id="eggnog_protein",
            ... )
            >>> selection.extract_mapping().select(
            ...     "GroupId", "CogCategory"
            ... ).rows()
            [('up', 'E'), ('up', 'G'), ('up', 'S')]
        """
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
        """Reject an unstable lazy build for the SQLite-backed dataset.

        eggNOG mapping expands a SQLite table into a flat TSV before Polars can
        scan it lazily, so callers should use `write_tidy()` for full writes.

        Raises:
            NotImplementedError: Always. The intermediate TSV cannot outlive a
                standalone lazy dataset handle.

        Examples:
            Use `write_tidy()` instead of requesting an unstable lazy dataset:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> db.build_tidy()
            Traceback (most recent call last):
            ...
            NotImplementedError: EggnogDb.build_tidy() cannot return a stable lazy dataset without a materialized mapping source; use EggnogDb.write_tidy().
            >>> [asset.path for asset in db.write_tidy("out/eggnog").assets]
            ['mapping.parquet']
        """
        raise NotImplementedError(
            "EggnogDb.build_tidy() cannot return a stable lazy dataset without "
            "a materialized mapping source; use EggnogDb.write_tidy()."
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the complete eggNOG mapping as canonical parquet.

        Args:
            dir_out: Destination directory for `mapping.parquet` and the
                optional manifest.
            should_write_manifest: Whether to write `manifest.json` and return
                its content in the report.
            should_hash_assets: Whether to calculate asset SHA-256 values for
                a requested manifest.

        Returns:
            A report describing the canonical mapping asset and optional
            manifest.

        Notes:
            The implementation streams SQLite rows through a temporary TSV.
            Wrapped databases and temporary assets use the configured
            `dir_tmp` when one was supplied at construction.

        Examples:
            Write the canonical mapping and inspect the published asset name:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> report = db.write_tidy("out/eggnog")
            >>> [asset.path for asset in report.assets]
            ['mapping.parquet']
        """
        with tempfile.TemporaryDirectory(
            prefix="bioextract-eggnog-",
            dir=None if self.snapshot.dir_tmp is None else self.snapshot.dir_tmp,
        ) as dir_tmp:
            file_mapping_tsv = Path(dir_tmp) / "mapping.tsv"
            write_mapping_tsv(
                file_eggnog_db=self.snapshot.file_eggnog_db,
                dir_tmp=self.snapshot.dir_tmp,
                df_cog_fun=self.read_cog_fun(),
                file_out=file_mapping_tsv,
            )
            dataset = EggnogTidyDataset(
                frames={"mapping": scan_mapping_tsv(file_mapping_tsv)},
                source=self._tidy_sources(),
                schema_version=SCHEMA_VERSION,
                build_id_prefix=f"eggnog-mapping-{self.snapshot.file_eggnog_db.stem}",
                assets=tuple(
                    TidyAsset(path=path, kind=kind, frame_name=frame_name)
                    for path, kind, frame_name in ASSET_SPECS
                ),
            )
            return dataset.write(
                Path(dir_out),
                should_write_manifest=should_write_manifest,
                should_hash_assets=should_hash_assets,
            )

    def read_cog_fun(self) -> pl.DataFrame:
        """Read and cache the optional COG function lookup.

        Returns:
            A frame keyed by `CogCategory`. If no lookup file was configured,
            the frame is empty but retains the expected lookup schema.

        Examples:
            Resolve a COG category to its display name:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db",
            ...     file_cog_fun="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.read_cog_fun().select(
            ...     "CogCategory", "CogName"
            ... ).head(1).rows()
            [('E', 'Amino acid transport and metabolism')]
        """
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
    """Materialize one single or grouped eggNOG mapping query.

    Instances are created by `EggnogDb.select_ids()` or
    `EggnogDb.select_groups()`. Mapping and unmapped frames are cached
    independently after first extraction.

    Examples:
        Inspect the COG categories retained by a selection:

        >>> db = EggnogDb.from_files(
        ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
        ... )
        >>> selection = db.select_ids(
        ...     ["9606.ENSP1"], kind_input_id="eggnog_protein"
        ... )
        >>> selection.extract_mapping().get_column("CogCategory").to_list()
        ['E', 'G', 'S']
    """

    dataset: EggnogDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: EggnogInputIdKind
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> db.select_ids(
            ...     ["9606.ENSP1"], kind_input_id="eggnog_protein"
            ... ).is_grouped
            False
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract mapping rows for the normalized selection.

        Returns:
            A frame prefixed by `InputId` and `KindInputId`; grouped selections
            additionally begin with `GroupId`. Many-to-many annotations remain
            expanded.

        Examples:
            Extract the normalized input alongside its COG categories:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(
            ...     ["9606.ENSP1"], kind_input_id="eggnog_protein"
            ... )
            >>> selection.extract_mapping().select(
            ...     "InputId", "CogCategory"
            ... ).rows()
            [('9606.ENSP1', 'E'), ('9606.ENSP1', 'G'), ('9606.ENSP1', 'S')]
        """
        if self._df_mapping is None:
            self._df_mapping = select_mapping_frame(
                file_eggnog_db=self.dataset.snapshot.file_eggnog_db,
                dir_tmp=self.dataset.snapshot.dir_tmp,
                df_input_ids=self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
                df_cog_fun=self.dataset.read_cog_fun(),
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no eggNOG mapping row.

        Returns:
            `InputId` for a single selection, or `GroupId, InputId` for a
            grouped selection.

        Examples:
            Report an identifier absent from the local snapshot:

            >>> db = EggnogDb.from_files(
            ...     file_eggnog_db="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(
            ...     ["9606.MISSING"], kind_input_id="eggnog_protein"
            ... )
            >>> selection.extract_unmapped_input_ids().to_dicts()
            [{'InputId': '9606.MISSING'}]
        """
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
