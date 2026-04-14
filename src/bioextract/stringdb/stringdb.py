import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
    validate_required_cols,
)
from bioextract.stringdb.spec import StringResourceLimits

from .constant import (
    DEFAULT_SOURCE_RANK_MAP,
    SCHEMA_ALIASES,
    SCHEMA_EDGES,
    SCHEMA_GROUP_EDGES,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUP_STRING_MAPPING,
    SCHEMA_GROUPS,
    SCHEMA_LINKS,
    SCHEMA_PROTEIN_MAP,
    SCHEMA_UNMAPPED,
    StringDbVersion,
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
    "StringDb",
]


@dataclass(frozen=True, slots=True)
class _StringSnapshot:
    version: StringDbVersion
    file_aliases: Path | None = None
    file_links: Path | None = None


@dataclass(frozen=True, slots=True)
class _StringAliasInfo:
    col_string_id: str
    has_source: bool
    file_alias: Path


@dataclass(slots=True)
class StringDb:
    """Path-first access to a local STRING snapshot.

    `StringDb` is the public entrypoint for extracting STRING mappings and
    interaction edges from local `protein.aliases` and `protein.links` files.
    It keeps dataset-level configuration such as source-rank policy, supported
    STRING version, and fail-fast resource limits.

    Construct instances with :meth:`from_files` and then choose one of two
    query styles:

    - single selection via :meth:`select_ids`
    - grouped selection via :meth:`select_groups`

    Both query styles preserve lazy scan behavior for the large STRING inputs
    until the final materialized result is requested.

    Examples:
        Single query:

            db = StringDb.from_files(
                file_aliases="9606.protein.aliases.v12.0.txt.gz",
                file_links="9606.protein.links.v12.0.txt.gz",
            )
            df_edges = (
                db.select_ids(["TP53", "EGFR"])
                .with_score_min(400)
                .extract_edges()
            )

        Grouped query:

            db = StringDb.from_files(
                file_aliases="9606.protein.aliases.v12.0.txt.gz",
                file_links="9606.protein.links.v12.0.txt.gz",
            )
            df_edges = (
                db.select_groups({"TumorA": ["TP53", "EGFR"], "TumorB": ["CDK2"]})
                .with_score_min(400)
                .extract_edges()
            )
    """

    snapshot: _StringSnapshot
    source_rank_map: Mapping[str, int]
    limits: StringResourceLimits = field(default_factory=StringResourceLimits)
    _alias_schema_cached: _StringAliasInfo | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_SOURCE_RANK_MAP = DEFAULT_SOURCE_RANK_MAP
    DEFAULT_RESOURCE_LIMITS = StringResourceLimits()

    @classmethod
    def from_files(
        cls,
        file_aliases: os.PathLike[str] | None = None,
        file_links: os.PathLike[str] | None = None,
        *,
        source_rank_map: Mapping[str, int] = DEFAULT_SOURCE_RANK_MAP,
        limits: StringResourceLimits | None = None,
        version: StringDbVersion = "v12.0",
    ) -> StringDb:
        """Create a dataset handle from local STRING aliases and links files.

        Args:
            file_aliases: Path to a local STRING `protein.aliases` text or gzip
                file.
            file_links: Path to a local STRING `protein.links` text or gzip
                file.
            source_rank_map: Source-priority mapping used to break ties when the
                same input ID maps to the same STRING ID through multiple alias
                sources.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used.
            version: Declared STRING schema version used for CSV schema
                overrides and column validation.

        Returns:
            A dataset handle that can produce single or grouped selections.

        Raises:
            FileNotFoundError: If either input file does not exist.
            ValueError: If the requested version is unsupported or a configured
                file-size limit is exceeded.
        """
        if version not in SCHEMA_ALIASES or version not in SCHEMA_LINKS:
            raise ValueError(f"Unsupported STRING version: {version}")

        limits_resolved = StringResourceLimits() if limits is None else limits

        if file_aliases is not None:
            file_aliases = Path(file_aliases)
            if not file_aliases.exists():
                raise FileNotFoundError(
                    f"STRING aliases file not found: {file_aliases}"
                )
            validate_file_size(
                file_path=file_aliases,
                size_max=limits_resolved.file_aliases_bytes_max,
                label="STRING aliases file",
            )

        if file_links is not None:
            file_links = Path(file_links)
            if not file_links.exists():
                raise FileNotFoundError(f"STRING links file not found: {file_links}")
            validate_file_size(
                file_path=file_links,
                size_max=limits_resolved.file_links_bytes_max,
                label="STRING links file",
            )

        return cls(
            snapshot=_StringSnapshot(
                version=version,
                file_aliases=file_aliases,
                file_links=file_links,
            ),
            source_rank_map=dict(source_rank_map),
            limits=limits_resolved,
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

        Raises:
            ValueError: If the normalized input-ID count exceeds the configured
                dataset limits.
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return StringSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            thr_score_min=0,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
    ) -> StringSelection:
        """Create a grouped selection from multiple input-ID sets.

        Each group key is normalized with `str(...).strip()`. Input IDs within
        each group follow the same normalization rules as :meth:`select_ids`.
        Grouped extraction keeps groups isolated in the returned flat tables by
        carrying a `GroupId` column through mapping, unmapped, edge, and metric
        outputs.

        Args:
            group_to_ids: Mapping from group label to iterable of input
                identifiers.

        Returns:
            A `StringSelection` in grouped mode. The returned selection emits
            flat grouped tables with `GroupId` as the leading column.

        Raises:
            ValueError: If any normalized `GroupId` is empty, if normalized
                `GroupId` values are duplicated, if the number of groups
                exceeds the configured limit, or if the total normalized
                input-ID count exceeds the configured limit.
        """
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
        return StringSelection(
            dataset=self,
            _df_groups=grp_in_frames.df_groups,
            _df_input_ids=grp_in_frames.df_input_ids,
            thr_score_min=0,
        )

    @property
    def alias_schema(self) -> _StringAliasInfo | None:
        if self._alias_schema_cached is not None:
            return self._alias_schema_cached

        if self.snapshot.file_aliases is None:
            return None

        lf_aliases = scan_aliases(
            self.snapshot.file_aliases, version=self.snapshot.version
        )
        cols_available = lf_aliases.collect_schema().names()
        col_id = infer_alias_id_col(
            cols_available,
            version=self.snapshot.version,
        )
        validate_alias_required_cols(cols_available, version=self.snapshot.version)
        self._alias_schema_cached = _StringAliasInfo(
            col_string_id=col_id,
            has_source="source" in cols_available,
            file_alias=self.snapshot.file_aliases,
        )

        return self._alias_schema_cached


@dataclass(slots=True)
class StringSelection:
    """Selection handle for both single and grouped STRING queries.

    `StringSelection` is the only public selection type exposed by
    `bioextract.stringdb`. Its output schemas depend on mode:

    - selections created by :meth:`StringDb.select_ids` return single-query
      tables without `GroupId`
    - selections created by :meth:`StringDb.select_groups` return grouped flat
      tables with leading `GroupId`
    """

    dataset: StringDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    thr_score_min: int = 0
    _df_protein_map: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)
    _df_string_ids: pl.DataFrame | None = field(default=None, repr=False)
    _df_edges: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Returns:
            `True` when the selection was created by :meth:`StringDb.select_groups`;
            otherwise `False`.
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

    def with_score_min(self, thr_score_min: int) -> StringSelection:
        """Create a new selection with a different minimum STRING score.

        Cached mapping-related frames are reused. Edge-related caches are not
        reused because the score threshold changes the edge result.

        Args:
            thr_score_min: Minimum `combined_score` required for retained
                STRING edges.

        Returns:
            A new selection sharing cached mapping state with the current
            selection.
        """
        return StringSelection(
            dataset=self.dataset,
            _df_input_ids=self._df_input_ids,
            _df_groups=self._df_groups,
            thr_score_min=int(thr_score_min),
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
        """
        if self.dataset.alias_schema is None:
            raise ValueError("Cannot extract STRING mapping without aliases file")

        if self._df_protein_map is not None:
            return self._df_protein_map

        if self._df_input_ids.height == 0:
            self._df_protein_map = pl.DataFrame(schema=self._schema_string_map_result)
            return self._df_protein_map

        lf_aliases = scan_aliases(
            self.dataset.alias_schema.file_alias, version=self.dataset.snapshot.version
        )
        lf_input_ids = self._df_input_ids.lazy()
        col_group_id = list(self._col_group_id)
        self._df_protein_map = (
            create_string_mapping_lazy_frame(
                lf_aliases=lf_aliases,
                lf_input_ids=lf_input_ids,
                source_rank_map=self.dataset.source_rank_map,
                col_string_id_aliases=self.dataset.alias_schema.col_string_id,
                has_source_aliases=self.dataset.alias_schema.has_source,
                cols_partition=col_group_id + ["InputId", "StringId"],
                cols_sort_prefix=col_group_id + ["InputId", "StringId"],
                cols_select_out=col_group_id + ["InputId", "StringId", "MapSource"],
            )
            .sort(col_group_id + ["InputId", "StringId", "MapSource"])
            .collect()
        )

        return self._df_protein_map

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that were not mapped to STRING IDs.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `InputId`
            - grouped selection: `GroupId`, `InputId`

            Each row represents a normalized input ID that did not resolve
            through the aliases table.
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

        version_link = self.dataset.snapshot.version
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
                thr_score_min=self.thr_score_min,
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
