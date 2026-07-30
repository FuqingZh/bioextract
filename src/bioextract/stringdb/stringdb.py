from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_required_cols,
)

from .constant import (
    DEFAULT_SOURCE_RANK_MAP,
    SCHEMA_EDGES,
    SCHEMA_GROUP_EDGES,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUP_STRING_MAPPING,
    SCHEMA_GROUPS,
    SCHEMA_LINKS,
    SCHEMA_PROTEIN_MAP,
    SCHEMA_UNMAPPED,
    StringDatabaseVersion,
)
from .util import (
    create_edges_lazy_frame,
    create_string_mapping_lazy_frame,
    infer_alias_id_col,
    scan_aliases,
    scan_links,
    validate_alias_required_cols,
)

__all__ = [
    "STRINGDatabase",
]


@dataclass(frozen=True, slots=True)
class _StringSnapshot:
    parser_version: StringDatabaseVersion
    release_version: str | None = None
    file_aliases: Path | None = None
    file_links: Path | None = None


@dataclass(frozen=True, slots=True)
class _StringAliasInfo:
    col_string_id: str
    has_source: bool
    file_alias: Path


@dataclass(slots=True)
class STRINGDatabase:
    """Path-first access to a local STRING snapshot.

    `STRINGDatabase` is the public entrypoint for extracting STRING mappings and
    interaction edges from local `protein.aliases` and `protein.links` files.
    It keeps dataset-level configuration such as source-rank policy and the
    supported STRING version.

    Construct instances with :meth:`from_files` and then choose one of two
    query styles:

    - single selection via :meth:`select_ids`
    - grouped selection via :meth:`select_groups`

    Both query styles preserve lazy scan behavior for the large STRING inputs
    until the final materialized result is requested.

    Examples:
        Create one snapshot handle and run a single query:

        >>> db = STRINGDatabase.from_files(
        ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz",
        ...     links="fixtures/string/9606.protein.links.v12.0.txt.gz",
        ... )
        >>> (
        ...     db.select_ids(["TP53", "EGFR"])
        ...     .with_min_combined_score(400)
        ...     .extract_edges()
        ...     .to_dicts()
        ... )
        [{'StringIdA': '9606.ENSP0001', 'StringIdB': '9606.ENSP0002', 'Score': 700}]
    """

    snapshot: _StringSnapshot
    source_rank_map: Mapping[str, int]
    _alias_schema_cached: _StringAliasInfo | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_SOURCE_RANK_MAP = DEFAULT_SOURCE_RANK_MAP

    @classmethod
    def from_files(
        cls,
        aliases: os.PathLike[str] | None = None,
        links: os.PathLike[str] | None = None,
        *,
        rank_by_source: Mapping[str, int] = DEFAULT_SOURCE_RANK_MAP,
        release_version: str | None = None,
    ) -> STRINGDatabase:
        """Create a dataset handle from local STRING aliases and links files.

        Args:
            aliases: Path to a local STRING `protein.aliases` text or gzip
                file.
            links: Path to a local STRING `protein.links` text or gzip
                file.
            rank_by_source: Source-priority mapping used to break ties when the
                same input ID maps to the same STRING ID through multiple alias
                sources.
            release_version: Optional official STRING release identity supplied
                by the caller. The parser profile is selected internally and
                every input header is validated against it.

        Returns:
            A dataset handle that can produce single or grouped selections.

        Raises:
            FileNotFoundError: If either input file does not exist.
            ValueError: If the release identity is empty or an input does not
                match the supported content profile.

        Examples:
            Open a STRING fixture snapshot:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz",
            ...     links="fixtures/string/9606.protein.links.v12.0.txt.gz",
            ... )
            >>> (
            ...     db.select_ids(["TP53"])
            ...     .extract_string_mapping()
            ...     .select("InputId", "StringId")
            ...     .to_dicts()
            ... )
            [{'InputId': 'TP53', 'StringId': '9606.ENSP0001'}]
        """
        if release_version is not None and not str(release_version).strip():
            raise ValueError("STRING release_version must be non-empty")

        file_aliases = aliases
        file_links = links
        if file_aliases is not None:
            file_aliases = Path(file_aliases)
            if not file_aliases.exists():
                raise FileNotFoundError(
                    f"STRING aliases file not found: {file_aliases}"
                )

        if file_links is not None:
            file_links = Path(file_links)
            if not file_links.exists():
                raise FileNotFoundError(f"STRING links file not found: {file_links}")

        return cls(
            snapshot=_StringSnapshot(
                parser_version="v12.0",
                release_version=release_version,
                file_aliases=file_aliases,
                file_links=file_links,
            ),
            source_rank_map=dict(rank_by_source),
        )

    def select_ids(self, ids: Iterable[str]) -> StringSelection:
        """Create a single-query selection from input IDs.

        Input IDs are normalized before selection:

        - surrounding whitespace is stripped
        - UniProt pipe-style values such as `sp|P04637|P53_HUMAN` are reduced
          to the middle accession token
        - empty normalized IDs are dropped
        - duplicates are removed

        Args:
            ids: Input protein, gene, or alias identifiers to resolve
                against the STRING aliases table.

        Returns:
            A `StringSelection` in single-query mode. The returned selection
            can extract mapping, unmapped IDs, and edges.

        Examples:
            Normalize an input ID and retain an unmapped alias:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> selection = db.select_ids([" TP53 ", "MISSING"])
            >>> selection.extract_string_mapping()["InputId"].to_list()
            ['TP53']
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return StringSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            min_combined_score=0,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> StringSelection:
        """Create a grouped selection from multiple input-ID sets.

        Each group key is normalized with `str(...).strip()`. Input IDs within
        each group follow the same normalization rules as :meth:`select_ids`.
        Grouped extraction keeps groups isolated in the returned flat tables by
        carrying a `GroupId` column through mapping, unmapped, edge, and metric
        outputs.

        Args:
            ids_by_group: Mapping from group label to iterable of input
                identifiers.

        Returns:
            A `StringSelection` in grouped mode. The returned selection emits
            flat grouped tables with `GroupId` as the leading column.

        Raises:
            ValueError: If any normalized `GroupId` is empty, if normalized
                `GroupId` values are duplicated.

        Examples:
            Keep each mapping in its original comparison:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> (
            ...     db.select_groups({"up": ["TP53"], "down": ["EGFR"]})
            ...     .extract_string_mapping()
            ...     .select("GroupId", "InputId")
            ...     .to_dicts()
            ... )
            [{'GroupId': 'down', 'InputId': 'EGFR'}, {'GroupId': 'up', 'InputId': 'TP53'}]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return StringSelection(
            dataset=self,
            _df_groups=grp_in_frames.df_groups,
            _df_input_ids=grp_in_frames.df_input_ids,
            min_combined_score=0,
        )

    @property
    def _alias_schema(self) -> _StringAliasInfo | None:
        """Return and cache aliases-file parsing metadata.

        Returns:
            Schema information required by mapping extraction, or ``None``
            when the snapshot has no aliases file.

        Raises:
            ValueError: If an aliases file is present but its required columns
                do not match the declared STRING version.

        Notes:
            Mapping extraction owns the public output schema. Keep this
            version-specific parser metadata private because its return type
            and source-column names are implementation details.
        """
        if self._alias_schema_cached is not None:
            return self._alias_schema_cached

        if self.snapshot.file_aliases is None:
            return None

        lf_aliases = scan_aliases(
            self.snapshot.file_aliases,
            version=self.snapshot.parser_version,
        )
        cols_available = lf_aliases.collect_schema().names()
        col_id = infer_alias_id_col(
            cols_available,
            version=self.snapshot.parser_version,
        )
        validate_alias_required_cols(
            cols_available, version=self.snapshot.parser_version
        )
        self._alias_schema_cached = _StringAliasInfo(
            col_string_id=col_id,
            has_source="source" in cols_available,
            file_alias=self.snapshot.file_aliases,
        )

        return self._alias_schema_cached


@dataclass(slots=True)
class StringSelection:
    """Selection handle for both single and grouped STRING queries.

    `StringSelection` is returned by :meth:`STRINGDatabase.select_ids` and
    :meth:`STRINGDatabase.select_groups`. Its output schemas depend on mode:

    - selections created by :meth:`STRINGDatabase.select_ids` return single-query
      tables without `GroupId`
    - selections created by :meth:`STRINGDatabase.select_groups` return grouped flat
      tables with leading `GroupId`

    Examples:
        Use a returned selection to resolve an input alias:

        >>> db = STRINGDatabase.from_files(
        ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
        ... )
        >>> selection = db.select_ids(["TP53"])
        >>> (
        ...     selection.extract_string_mapping()
        ...     .select("InputId", "StringId")
        ...     .to_dicts()
        ... )
        [{'InputId': 'TP53', 'StringId': '9606.ENSP0001'}]
    """

    dataset: STRINGDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    min_combined_score: int = 0
    _df_protein_map: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)
    _df_string_ids: pl.DataFrame | None = field(default=None, repr=False)
    _df_edges: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Returns:
            `True` when the selection was created by :meth:`STRINGDatabase.select_groups`;
            otherwise `False`.

        Examples:
            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> db.select_groups({"up": ["TP53"]}).is_grouped
            True
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the name of the group ID column when in grouped mode.

        Returns:
            A singleton tuple containing the name of the group ID column when
            `is_grouped` is `True`; otherwise an empty tuple.
        """
        return ("GroupId",) if self.is_grouped else ()

    @property
    def _schema_string_map_result(self) -> SchemaDict:
        """Return the expected schema for the mapping table output."""
        return SCHEMA_GROUP_STRING_MAPPING if self.is_grouped else SCHEMA_PROTEIN_MAP

    @property
    def _schema_edges_result(self) -> SchemaDict:
        """Return the expected schema for the edges table output."""
        return SCHEMA_GROUP_EDGES if self.is_grouped else SCHEMA_EDGES

    def with_min_combined_score(self, min_combined_score: int) -> StringSelection:
        """Create a new selection with a different minimum STRING score.

        Cached mapping-related frames are reused. Edge-related caches are not
        reused because the score threshold changes the edge result.

        Args:
            min_combined_score: Minimum `combined_score` required for retained
                STRING edges.

        Returns:
            A new selection sharing cached mapping state with the current
            selection.

        Examples:
            Remove an edge whose combined score is below 700:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz",
            ...     links="fixtures/string/9606.protein.links.v12.0.txt.gz",
            ... )
            >>> selection = db.select_ids(["TP53", "EGFR", "CDK2"])
            >>> (
            ...     selection.with_min_combined_score(400)
            ...     .extract_edges()["Score"]
            ...     .sort()
            ...     .to_list()
            ... )
            [450, 700]
            >>> selection.with_min_combined_score(700).extract_edges()["Score"].to_list()
            [700]
        """
        return StringSelection(
            dataset=self.dataset,
            _df_input_ids=self._df_input_ids,
            _df_groups=self._df_groups,
            min_combined_score=int(min_combined_score),
            _df_protein_map=self._df_protein_map,
            _df_unmapped=self._df_unmapped,
            _df_string_ids=self._df_string_ids,
        )

    def extract_string_mapping(self) -> pl.DataFrame:
        """Extract the input-to-STRING mapping table for this selection.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `InputId`, `StringId`, `MapSource`
            - grouped selection: `GroupId`, `InputId`, `StringId`, `MapSource`

        Raises:
            ValueError: If the aliases file is missing required columns for the
                configured STRING version.

        Examples:
            Resolve a gene name and report the chosen alias source:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> db.select_ids(["TP53"]).extract_string_mapping().to_dicts()
            [{'InputId': 'TP53', 'StringId': '9606.ENSP0001', 'MapSource': 'UniProt_GN_Name'}]
        """
        if self.dataset._alias_schema is None:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            raise ValueError("Cannot extract STRING mapping without aliases file")

        if self._df_protein_map is not None:
            return self._df_protein_map

        if self._df_input_ids.height == 0:
            self._df_protein_map = pl.DataFrame(schema=self._schema_string_map_result)
            return self._df_protein_map

        lf_aliases = scan_aliases(
            self.dataset._alias_schema.file_alias,  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            version=self.dataset.snapshot.parser_version,
        )
        lf_input_ids = self._df_input_ids.lazy()
        col_group_id = list(self._col_group_id)
        self._df_protein_map = (
            create_string_mapping_lazy_frame(
                lf_aliases=lf_aliases,
                lf_input_ids=lf_input_ids,
                source_rank_map=self.dataset.source_rank_map,
                col_string_id_aliases=self.dataset._alias_schema.col_string_id,  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                has_source_aliases=self.dataset._alias_schema.has_source,  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                cols_partition=col_group_id + ["InputId", "StringId"],
                cols_sort_prefix=col_group_id + ["InputId", "StringId"],
                cols_select_out=col_group_id + ["InputId", "StringId", "MapSource"],
            )
            .sort(col_group_id + ["InputId", "StringId", "MapSource"])
            .collect()
        )

        return self._df_protein_map

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that were not mapped to STRING IDs.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `InputId`
            - grouped selection: `GroupId`, `InputId`

            Each row represents a normalized input ID that did not resolve
            through the aliases table.

        Examples:
            Report an identifier absent from the aliases fixture:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> db.select_ids(["MISSING"]).extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            cols_index = list(self._col_group_id) + ["InputId"]
            df_mapped_input_ids = self.extract_string_mapping().select(cols_index)
            self._df_unmapped = (
                self._df_input_ids.join(
                    df_mapped_input_ids.unique(),
                    on=cols_index,
                    how="anti",
                )
                .select(cols_index)
                .sort(cols_index)
            )

        return self._df_unmapped

    def extract_edges(self) -> pl.DataFrame:
        """Extract the STRING subnetwork induced by the mapped STRING IDs.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `StringIdA`, `StringIdB`, `Score`
            - grouped selection: `GroupId`, `StringIdA`, `StringIdB`, `Score`

        Raises:
            ValueError: If the links file is missing required columns for the
                configured STRING version.

        Examples:
            Extract the induced edge and its combined score:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz",
            ...     links="fixtures/string/9606.protein.links.v12.0.txt.gz",
            ... )
            >>> db.select_ids(["TP53", "EGFR"]).extract_edges().to_dicts()
            [{'StringIdA': '9606.ENSP0001', 'StringIdB': '9606.ENSP0002', 'Score': 700}]
        """
        if self.dataset.snapshot.file_links is None:
            raise ValueError("Cannot extract STRING edges without links file")
        if self._df_edges is not None:
            return self._df_edges

        col_group_id = list(self._col_group_id)
        df_string_ids = self._extract_string_ids()
        if df_string_ids.height == 0:
            self._df_edges = pl.DataFrame(schema=self._schema_edges_result)
            return self._df_edges
        lf_string_ids = df_string_ids.lazy()

        version_link = self.dataset.snapshot.parser_version
        lf_links = scan_links(self.dataset.snapshot.file_links, version=version_link)
        validate_required_cols(
            cols_available=lf_links.collect_schema().names(),
            cols_required=SCHEMA_LINKS[version_link].keys(),
            context=f"STRING links file for version {version_link}",
        )
        self._df_edges = (
            create_edges_lazy_frame(
                lf_links=lf_links,
                lf_string_ids_a=lf_string_ids.rename({"StringId": "StringIdA"}),
                lf_string_ids_b=lf_string_ids.rename({"StringId": "StringIdB"}),
                min_combined_score=self.min_combined_score,
                cols_join_left_a="StringIdA",
                cols_join_right_a="StringIdA",
                cols_join_left_b=col_group_id + ["StringIdB"],
                cols_join_right_b=col_group_id + ["StringIdB"],
                cols_partition=col_group_id + ["_Lo", "_Hi"],
                cols_select_out=col_group_id + ["StringIdA", "StringIdB", "Score"],
            )
            .sort(col_group_id + ["StringIdA", "StringIdB"])
            .collect()
        )

        return self._df_edges

    def _extract_string_ids(self) -> pl.DataFrame:
        cols_select = list(self._col_group_id) + ["StringId"]
        if self._df_string_ids is None:
            self._df_string_ids = (
                self.extract_string_mapping()
                .select(cols_select)
                .unique()
                .sort(cols_select)
            )
        return self._df_string_ids
