from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._publication import ParquetWriteResult
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource

from .constant import (
    ASSET_SPECS,
    MEDIA_TYPE_SQLITE,
    MEDIA_TYPE_SQLITE_GZIP,
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
    EggnogNamespace,
)
from .util import (
    build_mapping_frame,
    extract_unmatched_ids_frame,
    read_cog_fun_frame,
    scan_mapping_tsv,
    select_mapping_frame,
    write_mapping_tsv,
)

__all__ = [
    "EggNOGDatabase",
    "EggnogSelection",
    "EggnogTidyDataset",
]


@dataclass(frozen=True, slots=True)
class _EggnogSnapshot:
    file_eggnog_db: Path
    file_cog_fun: Path | None
    dir_tmp: Path | None


EggnogTidyDataset = TidyDataset


@dataclass(slots=True)
class EggNOGDatabase:
    """Access one local eggNOG mapper resource snapshot.

    Construction validates paths without expanding the SQLite mapping. The
    optional COG function table enriches `CogClass` and
    `CogName`; those columns remain null when the table is omitted. Materialized
    mapping and lookup frames are cached on the handle.

    Examples:
        Read one enriched protein-to-COG annotation:

        >>> db = EggNOGDatabase.from_files(
        ...     eggnog_database="fixtures/eggnog/eggnog.db",
        ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
        ... )
        >>> db.extract_mapping().select(
        ...     "EggnogProteinId", "CogCategory", "CogName"
        ... ).head(1).rows()
        [('9606.ENSP1', 'E', 'Amino acid transport and metabolism')]
    """

    snapshot: _EggnogSnapshot
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_cog_fun: pl.DataFrame | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_files(
        cls,
        *,
        eggnog_database: os.PathLike[str] | str,
        cog_functions: os.PathLike[str] | str | None = None,
        temp_dir: os.PathLike[str] | str | None = None,
    ) -> EggNOGDatabase:
        """Create a handle from explicit eggNOG mapper files.

        Args:
            eggnog_database: Path to an eggNOG mapper SQLite database or its
                `.gz` wrapper.
            cog_functions: Optional COG function table used to populate class
                and display-name columns.
            temp_dir: Optional scratch directory for decompressing a wrapped
                SQLite database and staging full tidy exports.

        Returns:
            A lightweight handle that defers database reads until extraction.

        Raises:
            FileNotFoundError: If a supplied source file does not exist.

        Examples:
            Enrich COG categories with names from the optional lookup:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db",
            ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.extract_mapping().select(
            ...     "CogCategory", "CogName"
            ... ).head(1).rows()
            [('E', 'Amino acid transport and metabolism')]
        """
        file_eggnog_db = _validate_file(
            eggnog_database,
            label="eggNOG SQLite database file",
        )
        file_cog_fun = cog_functions
        if file_cog_fun is not None:
            file_cog_fun = _validate_file(
                file_cog_fun,
                label="COG function table file",
            )
        dir_tmp_resolved = Path(temp_dir) if temp_dir is not None else None
        if dir_tmp_resolved is not None:
            dir_tmp_resolved.mkdir(parents=True, exist_ok=True)

        return cls(
            snapshot=_EggnogSnapshot(
                file_eggnog_db=file_eggnog_db,
                file_cog_fun=file_cog_fun,
                dir_tmp=dir_tmp_resolved,
            ),
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full protein-to-COG mapping table.

        Returns:
            A cached frame that preserves protein-to-OG and OG-to-COG
            many-to-many relationships. Missing optional COG metadata remains
            null.

        Examples:
            Preserve one row for each protein, OG, and COG-category match:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db",
            ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
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
    ) -> EggnogSelection:
        """Create a selection for one normalized ID set.

        Args:
            ids: eggNOG protein IDs. Empty values are discarded and duplicate
                normalized IDs collapse to one input row.
        Returns:
            A selection handle whose outputs include `InputId` and
            `InputNamespace` provenance columns.

        Examples:
            Return annotations only for the selected protein ID:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.ENSP1"])
            >>> selection.extract_mapping().select(
            ...     "InputId", "EggnogOgId", "CogCategory"
            ... ).rows()
            [('9606.ENSP1', 'OG0001', 'E'), ('9606.ENSP1', 'OG0001', 'G'), ('9606.ENSP1', 'OG0002', 'S')]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return EggnogSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> EggnogSelection:
        """Create a selection that preserves caller-defined groups.

        Args:
            ids_by_group: Mapping from unique, non-empty group IDs to eggNOG
                protein IDs.
        Returns:
            A selection handle whose mapping and unmapped outputs retain
            `GroupId`.

        Raises:
            ValueError: If group IDs are invalid.

        Examples:
            Preserve a comparison label on every selected annotation row:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["9606.ENSP1"], "down": ["9606.MISSING"]},
            ... )
            >>> selection.extract_mapping().select(
            ...     "GroupId", "CogCategory"
            ... ).rows()
            [('up', 'E'), ('up', 'G'), ('up', 'S')]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return EggnogSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
        )

    def write_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Stream the complete eggNOG mapping into one atomic Parquet file.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_parquet(Path(dir_out) / "eggnog.parquet")
            ...     result.resource_name.startswith("eggnog-")
            True
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
                path=file_mapping_tsv,
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
                resource_name="eggnog",
            )
            return dataset.write_parquet(path, if_exists=if_exists)

    def read_cog_fun(self) -> pl.DataFrame:
        """Read and cache the optional COG function lookup.

        Returns:
            A frame keyed by `CogCategory`. If no lookup file was configured,
            the frame is empty but retains the expected lookup schema.

        Examples:
            Resolve a COG category to its display name:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db",
            ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
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

    Instances are created by `EggNOGDatabase.select_ids()` or
    `EggNOGDatabase.select_groups()`. Mapping and unmapped frames are cached
    independently after first extraction.

    Examples:
        Inspect the COG categories retained by a selection:

        >>> db = EggNOGDatabase.from_files(
        ...     eggnog_database="fixtures/eggnog/eggnog.db"
        ... )
        >>> selection = db.select_ids(["9606.ENSP1"])
        >>> selection.extract_mapping().get_column("CogCategory").to_list()
        ['E', 'G', 'S']
    """

    dataset: EggNOGDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    namespace: EggnogNamespace = field(
        default="eggnog_protein",
        init=False,
    )
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> db.select_ids(["9606.ENSP1"]).is_grouped
            False
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract mapping rows for the normalized selection.

        Returns:
            A frame prefixed by `InputId` and `InputNamespace`; grouped selections
            additionally begin with `GroupId`. Many-to-many annotations remain
            expanded.

        Examples:
            Extract the normalized input alongside its COG categories:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.ENSP1"])
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
                namespace=self.namespace,
                cols_group_id=self._col_group_id,
                df_cog_fun=self.dataset.read_cog_fun(),
            )
        return self._df_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no eggNOG mapping row.

        Returns:
            `InputId` for a single selection, or `GroupId, InputId` for a
            grouped selection.

        Examples:
            Report an identifier absent from the local snapshot:

            >>> db = EggNOGDatabase.from_files(
            ...     eggnog_database="fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': '9606.MISSING'}]
        """
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                cols_group_id=self._col_group_id,
            )
        return self._df_unmapped


def _validate_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    return file_path
