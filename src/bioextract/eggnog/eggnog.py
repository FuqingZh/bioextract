from __future__ import annotations

import copy
import os
import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._lazy import register_deferred_frame_source
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)

from .constant import (
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    EggnogNamespace,
)
from .util import (
    extract_unmatched_ids_frame,
    is_gzip_file,
    read_cog_fun_frame,
    select_mapping_frame,
)

__all__ = [
    "EggNOGDatabase",
]


@dataclass(frozen=True, slots=True)
class _EggnogSnapshot:
    file_eggnog_db: Path
    file_cog_fun: Path | None
    dir_tmp: Path | None


@dataclass(slots=True)
class EggNOGDatabase:
    """Access one local eggNOG mapper resource snapshot.

    Construction validates paths without expanding the SQLite mapping. The
    optional COG function table enriches `cog_class` and
    `cog_name`; those columns remain null when the table is omitted. Materialized
    selected mapping frames are cached on their selection handles.

    Examples:
        Read one enriched protein-to-COG annotation:

        >>> db = EggNOGDatabase.from_sqlite(
        ...     "fixtures/eggnog/eggnog.db",
        ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
        ... )
        >>> db.select_ids(["9606.ENSP1"]).mappings().collect().select(
        ...     "name", "cog_category", "cog_name"
        ... ).head(1).rows()
        [('9606.ENSP1', 'E', 'Amino acid transport and metabolism')]
    """

    snapshot: _EggnogSnapshot
    _df_cog_fun: pl.DataFrame | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_sqlite(
        cls,
        source: os.PathLike[str] | str,
        *,
        cog_functions: os.PathLike[str] | str | None = None,
        temp_dir: os.PathLike[str] | str | None = None,
    ) -> EggNOGDatabase:
        """Create a handle from an eggNOG mapper SQLite database.

        Args:
            source: Path to an eggNOG mapper SQLite database or its gzip
                wrapper.
            cog_functions: Optional COG function table used to populate class
                and display-name columns.
            temp_dir: Optional scratch directory for temporarily decompressing
                a wrapped SQLite database.

        Returns:
            A lightweight handle that defers database reads until extraction.

        Raises:
            FileNotFoundError: If a supplied source file does not exist.

        Examples:
            Enrich COG categories with names from the optional lookup:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db",
            ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.select_ids(["9606.ENSP1"]).mappings().collect().select(
            ...     "cog_category", "cog_name"
            ... ).head(1).rows()
            [('E', 'Amino acid transport and metabolism')]
        """
        file_eggnog_db = _validate_file(
            source,
            label="eggNOG SQLite database file",
        )
        if is_gzip_file(file_eggnog_db):
            warnings.warn(
                "Compressed eggNOG SQLite source detected. bioextract must fully "
                "decompress it to a temporary SQLite file before access, and the "
                "decompressed file is not persisted. Repeated use may repeat this cost; "
                "for long-term use, decompress the source once and pass the .db file.",
                UserWarning,
                stacklevel=2,
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

    def select_ids(
        self,
        ids: Iterable[str],
    ) -> EggnogSelection:
        """Create a selection for one normalized ID set.

        Args:
            ids: eggNOG protein IDs. Empty values are discarded and duplicate
                normalized IDs collapse to one input row.
        Returns:
            A selection handle whose outputs include `input_id` and
            `input_namespace` provenance columns.

        Examples:
            Return annotations only for the selected protein ID:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.ENSP1"])
            >>> selection.mappings().select(
            ...     "input_id", "og", "cog_category"
            ... ).collect().rows()
            [('9606.ENSP1', 'OG0001', 'E'), ('9606.ENSP1', 'OG0001', 'G'), ('9606.ENSP1', 'OG0002', 'S')]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return EggnogSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
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
            `group_id`.

        Raises:
            ValueError: If group IDs are invalid.

        Examples:
            Preserve a comparison label on every selected annotation row:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["9606.ENSP1"], "down": ["9606.MISSING"]},
            ... )
            >>> selection.mappings().select(
            ...     "group_id", "cog_category"
            ... ).collect().rows()
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
            _df_group_membership=grp_in_frames.df_group_membership,
        )

    def cog_functions(self) -> pl.LazyFrame:
        """Return the optional COG function lookup lazily.

        Examples:
            >>> db.cog_functions().collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=dict.fromkeys(
                ["cog_category", "cog_class", "cog_name"],
                pl.String,
            ),
            frame=lambda: snapshot._read_cog_fun(),
        )

    def _read_cog_fun(self) -> pl.DataFrame:
        """Read and cache the optional COG function lookup.

        Returns:
            A frame keyed by `cog_category`. If no lookup file was configured,
            the frame is empty but retains the expected lookup schema.

        Examples:
            Resolve a COG category to its display name:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db",
            ...     cog_functions="fixtures/eggnog/cog-24.fun.tab",
            ... )
            >>> db.cog_functions().collect().select(
            ...     "cog_category", "cog_name"
            ... ).head(1).rows()
            [('E', 'Amino acid transport and metabolism')]
        """
        if self._df_cog_fun is None:
            self._df_cog_fun = read_cog_fun_frame(self.snapshot.file_cog_fun)
        return self._df_cog_fun


@dataclass(slots=True)
class EggnogSelection:
    """Materialize one single or grouped eggNOG mapping query.

    Instances are created by `EggNOGDatabase.select_ids()` or
    `EggNOGDatabase.select_groups()`. Mapping and unmapped frames are cached
    independently after first extraction.

    Examples:
        Inspect the COG categories retained by a selection:

        >>> db = EggNOGDatabase.from_sqlite(
        ...     "fixtures/eggnog/eggnog.db"
        ... )
        >>> selection = db.select_ids(["9606.ENSP1"])
        >>> selection.mappings().collect().get_column("cog_category").to_list()
        ['E', 'G', 'S']
    """

    dataset: EggNOGDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    namespace: EggnogNamespace = field(
        default="eggnog_protein",
        init=False,
    )
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `group_id` through outputs.

        Examples:
            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db"
            ... )
            >>> db.select_ids(["9606.ENSP1"]).is_grouped
            False
        """
        return self._df_groups is not None

    def mappings(self) -> pl.LazyFrame:
        """Return selected eggNOG mapping rows lazily.

        Examples:
            >>> selection.mappings().collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._df_group_membership is not None else []) + [
            "input_id",
            "input_namespace",
            "name",
            "og",
            "level",
            "description",
            "COG_categories",
            "cog_category",
            "cog_class",
            "cog_name",
        ]
        return register_deferred_frame_source(
            schema=dict.fromkeys(columns, pl.String),
            frame=snapshot._eager_mappings,
        )

    def _eager_mappings(self) -> pl.DataFrame:
        """Extract mapping rows for the normalized selection.

        Returns:
            A frame prefixed by `input_id` and `input_namespace`; grouped selections
            additionally begin with `group_id`. Many-to-many annotations remain
            expanded.

        Examples:
            Extract the normalized input alongside its COG categories:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.ENSP1"])
            >>> selection.mappings().select(
            ...     "input_id", "cog_category"
            ... ).collect().rows()
            [('9606.ENSP1', 'E'), ('9606.ENSP1', 'G'), ('9606.ENSP1', 'S')]
        """
        if self._df_mapping is None:
            self._df_mapping = select_mapping_frame(
                file_eggnog_db=self.dataset.snapshot.file_eggnog_db,
                dir_tmp=self.dataset.snapshot.dir_tmp,
                df_input_ids=self._df_input_ids,
                df_group_membership=self._df_group_membership,
                namespace=self.namespace,
                df_cog_fun=self.dataset._read_cog_fun(),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            )
        return self._df_mapping

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return selected IDs absent from the eggNOG mapping lazily.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (
            ["group_id", "input_id"]
            if self._df_group_membership is not None
            else ["input_id"]
        )
        return register_deferred_frame_source(
            schema=dict.fromkeys(columns, pl.String),
            frame=snapshot._eager_unmatched_ids,
        )

    def _eager_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no eggNOG mapping row.

        Returns:
            `input_id` for a single selection, or `group_id, input_id` for a
            grouped selection.

        Examples:
            Report an identifier absent from the local snapshot:

            >>> db = EggNOGDatabase.from_sqlite(
            ...     "fixtures/eggnog/eggnog.db"
            ... )
            >>> selection = db.select_ids(["9606.MISSING"])
            >>> selection.unmatched_ids().collect().to_dicts()
            [{'input_id': '9606.MISSING'}]
        """
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self._eager_mappings(),
                df_group_membership=self._df_group_membership,
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
