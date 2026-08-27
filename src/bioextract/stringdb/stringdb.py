from __future__ import annotations

import copy
import logging
import os
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import (
    _RelationScanRequest,  # pyright: ignore[reportPrivateUsage]  # typed source boundary
    register_replayable_source,
)
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    normalize_uniprot_selection_id,
    validate_required_cols,
)
from bioextract.errors import CapabilityError, IntegrityError

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
    StringNamespace,
)
from .util import (
    create_edges_lazy_frame,
    create_string_mapping_lazy_frame,
    infer_alias_id_col,
    scan_aliases,
    scan_links,
    validate_alias_required_cols,
)

logger = logging.getLogger(__name__)

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
        ...     .edges()
        ...     .collect().to_dicts()
        ... )
        [{'string_id_a': '9606.ENSP0001', 'string_id_b': '9606.ENSP0002', 'combined_score': 700}]
    """

    snapshot: _StringSnapshot
    source_rank_map: Mapping[str, int]
    _alias_schema_cached: _StringAliasInfo | None = field(
        default=None, init=False, repr=False
    )
    _taxa_validated: bool = field(default=False, init=False, repr=False)

    DEFAULT_SOURCE_RANK_MAP = DEFAULT_SOURCE_RANK_MAP

    @property
    def available_resources(self) -> tuple[str, ...]:
        """Return source roles available on this local snapshot.

        Examples:
            >>> database.available_resources  # doctest: +SKIP
            ('aliases', 'links')
        """
        resources: list[str] = []
        if self.snapshot.file_aliases is not None:
            resources.append("aliases")
        if self.snapshot.file_links is not None:
            resources.append("links")
        return tuple(resources)

    def scan_aliases(self) -> pl.LazyFrame:
        """Scan the official aliases relation in source column order.

        Raises:
            CapabilityError: If the required aliases source is absent from the
                local snapshot.

        Examples:
            >>> database.scan_aliases()  # doctest: +SKIP
            <LazyFrame ...>
        """
        if self.snapshot.file_aliases is None:
            raise CapabilityError("STRING aliases source is absent from this snapshot")
        return scan_aliases(
            self.snapshot.file_aliases,
            version=self.snapshot.parser_version,
        )

    def scan_links(self) -> pl.LazyFrame:
        """Scan the official links relation in source column order.

        Raises:
            CapabilityError: If the required links source is absent from the
                local snapshot.

        Examples:
            >>> database.scan_links()  # doctest: +SKIP
            <LazyFrame ...>
        """
        if self.snapshot.file_links is None:
            raise CapabilityError("STRING links source is absent from this snapshot")
        return scan_links(
            self.snapshot.file_links,
            version=self.snapshot.parser_version,
        )

    @classmethod
    def from_files(
        cls,
        *,
        aliases: os.PathLike[str] | str | None = None,
        links: os.PathLike[str] | str | None = None,
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
            ...     .mappings()
            ...     .select("input_id", "string_protein_id")
            ...     .collect().to_dicts()
            ... )
            [{'input_id': 'TP53', 'string_protein_id': '9606.ENSP0001', 'alias': 'TP53', 'source': 'UniProt_GN_Name'}]
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
            if file_aliases.suffix == ".gz":
                logger.warning(
                    "STRING aliases input %s is gzip-compressed; repeated scans "
                    "will re-decompress it. Uncompress long-lived snapshots for "
                    "better interactive performance.",
                    file_aliases,
                )

        if file_links is not None:
            file_links = Path(file_links)
            if not file_links.exists():
                raise FileNotFoundError(f"STRING links file not found: {file_links}")
            if file_links.suffix == ".gz":
                logger.warning(
                    "STRING links input %s is gzip-compressed; repeated scans "
                    "will re-decompress it. Uncompress long-lived snapshots for "
                    "better interactive performance.",
                    file_links,
                )

        return cls(
            snapshot=_StringSnapshot(
                parser_version="v12.0",
                release_version=release_version,
                file_aliases=file_aliases,
                file_links=file_links,
            ),
            source_rank_map=dict(rank_by_source),
        )

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: StringNamespace = "alias",
    ) -> StringSelection:
        """Create a single-query selection from input IDs.

        Input IDs are normalized before selection:

        - surrounding whitespace is stripped
        - alias inputs accept complete UniProt pipe-style values and reduce
          them to the accession
        - direct STRING namespace inputs retain pipe-bearing text exactly
        - empty normalized IDs are dropped
        - duplicates are removed

        Args:
            ids: Input protein, gene, or alias identifiers to resolve. With
                ``namespace='string'``, values are direct STRING protein IDs.

        Returns:
            A `StringSelection` in single-query mode. The returned selection
            can extract mapping, unmapped IDs, and edges.

        Raises:
            ValueError: If the namespace is unsupported or an alias
                pipe-bearing value is not one complete supported UniProt form.

        Examples:
            Normalize an input ID and retain an unmapped alias:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> selection = db.select_ids([" TP53 ", "MISSING"])
            >>> selection.mappings().select("input_id").collect().to_series().to_list()
            ['TP53']
            >>> selection.unmatched_ids().collect().to_dicts()
            [{'input_id': 'MISSING'}]
        """
        if namespace not in ("alias", "string"):
            raise ValueError(f"Unknown STRING namespace: {namespace!r}")
        normalized_ids = (
            (normalize_uniprot_selection_id(value) for value in ids)
            if namespace == "alias"
            else ids
        )
        df_input_ids = create_input_id_frame(
            normalized_ids, schema_unmapped=SCHEMA_UNMAPPED
        )
        return StringSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_group_membership=None,
            _df_groups=None,
            min_combined_score=0,
            namespace=namespace,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: StringNamespace = "alias",
    ) -> StringSelection:
        """Create a grouped selection from multiple input-ID sets.

        Each group key is normalized with `str(...).strip()`. Input IDs within
        each group follow the namespace-specific rules of :meth:`select_ids`.
        Alias resolution and source ranking run once per globally unique
        normalized ID, then the result is expanded through group membership.
        Edge extraction remains isolated by group. Returned mapping, unmatched,
        edge, and metric tables carry `group_id`.

        Args:
            ids_by_group: Mapping from group label to iterable of input
                identifiers.

        Returns:
            A `StringSelection` in grouped mode. The returned selection emits
            flat grouped tables with `group_id` as the leading column.

        Raises:
            ValueError: If any normalized `group_id` is empty, normalized
                `group_id` values are duplicated, or an alias pipe-bearing
                value is malformed.

        Examples:
            Keep each mapping in its original comparison:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
            ... )
            >>> (
            ...     db.select_groups({"up": ["TP53"], "down": ["EGFR"]})
            ...     .mappings()
            ...     .select("group_id", "input_id")
            ...     .collect().to_dicts()
            ... )
            [{'group_id': 'down', 'input_id': 'EGFR'}, {'group_id': 'up', 'input_id': 'TP53'}]
        """
        if namespace not in ("alias", "string"):
            raise ValueError(f"Unknown STRING namespace: {namespace!r}")
        normalized = (
            {
                group_id: (normalize_uniprot_selection_id(value) for value in input_ids)
                for group_id, input_ids in ids_by_group.items()
            }
            if namespace == "alias"
            else ids_by_group
        )
        grp_in_frames = create_group_input_frames(
            normalized,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return StringSelection(
            dataset=self,
            _df_groups=grp_in_frames.df_groups,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_group_membership=grp_in_frames.df_group_membership,
            min_combined_score=0,
            namespace=namespace,
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

    def _validate_taxon_compatibility(self) -> None:
        if (
            self._taxa_validated
            or self.snapshot.file_aliases is None
            or self.snapshot.file_links is None
        ):
            return
        alias_info = self._alias_schema
        if alias_info is None:
            return
        alias_taxa = set(
            self.scan_aliases()
            .select(
                pl.col(alias_info.col_string_id)
                .str.split(".")
                .list.first()
                .alias("taxon")
            )
            .drop_nulls()
            .collect()["taxon"]
            .to_list()
        )
        link_taxa = set(
            pl.concat(
                [
                    self.scan_links().select(pl.col("protein1").alias("protein")),
                    self.scan_links().select(pl.col("protein2").alias("protein")),
                ]
            )
            .select(pl.col("protein").str.split(".").list.first().alias("taxon"))
            .drop_nulls()
            .collect()["taxon"]
            .to_list()
        )
        if alias_taxa and link_taxa and alias_taxa.isdisjoint(link_taxa):
            raise IntegrityError(
                "STRING aliases and links resources have incompatible taxon "
                f"prefixes: aliases={sorted(alias_taxa)}, links={sorted(link_taxa)}"
            )
        self._taxa_validated = True


@dataclass(slots=True)
class StringSelection:
    """Selection handle for both single and grouped STRING queries.

    `StringSelection` is returned by :meth:`STRINGDatabase.select_ids` and
    :meth:`STRINGDatabase.select_groups`. Its output schemas depend on mode:

    - selections created by :meth:`STRINGDatabase.select_ids` return single-query
      tables without `group_id`
    - selections created by :meth:`STRINGDatabase.select_groups` return grouped flat
      tables with leading `group_id`

    Examples:
        Use a returned selection to resolve an input alias:

        >>> db = STRINGDatabase.from_files(
        ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz"
        ... )
        >>> selection = db.select_ids(["TP53"])
        >>> (
        ...     selection.mappings()
        ...     .select("input_id", "string_protein_id")
        ...     .collect().to_dicts()
        ... )
        [{'input_id': 'TP53', 'string_protein_id': '9606.ENSP0001', 'alias': 'TP53', 'source': 'UniProt_GN_Name'}]
    """

    dataset: STRINGDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    namespace: StringNamespace = "alias"
    min_combined_score: int = 0

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `group_id` through outputs.

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
        return ("group_id",) if self.is_grouped else ()

    @property
    def _schema_string_map_result(self) -> SchemaDict:
        """Return the expected schema for the mapping table output."""
        if self.dataset._alias_schema is None:  # pyright: ignore[reportPrivateUsage]
            return (
                SCHEMA_GROUP_STRING_MAPPING if self.is_grouped else SCHEMA_PROTEIN_MAP
            )
        columns = [
            "input_id",
            self.dataset._alias_schema.col_string_id,  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            "alias",
            "source",
        ]  # pyright: ignore[reportPrivateUsage]
        if self.is_grouped:
            columns = ["group_id", *columns]
        return dict.fromkeys(columns, pl.String)

    @property
    def _schema_edges_result(self) -> SchemaDict:
        """Return the expected schema for the edges table output."""
        return SCHEMA_GROUP_EDGES if self.is_grouped else SCHEMA_EDGES

    def with_min_combined_score(self, min_combined_score: int) -> StringSelection:
        """Create a new selection with a different minimum STRING score.

        The returned selection captures the immutable input and threshold;
        each relation remains a replayable lazy plan.

        Args:
            min_combined_score: Minimum `combined_score` required for retained
                STRING edges.

        Returns:
            A new selection with the same inputs and a new edge threshold.

        Examples:
            Remove an edge whose combined score is below 700:

            >>> db = STRINGDatabase.from_files(
            ...     aliases="fixtures/string/9606.protein.aliases.v12.0.txt.gz",
            ...     links="fixtures/string/9606.protein.links.v12.0.txt.gz",
            ... )
            >>> selection = db.select_ids(["TP53", "EGFR", "CDK2"])
            >>> (
            ...     selection.with_min_combined_score(400)
            ...     .edges().select("combined_score").collect()["combined_score"]
            ...     .sort().to_list()
            ... )
            [450, 700]
            >>> selection.with_min_combined_score(700).edges().select("combined_score").collect()["combined_score"].to_list()
            [700]
        """
        return StringSelection(
            dataset=self.dataset,
            _df_input_ids=self._df_input_ids,
            _df_group_membership=self._df_group_membership,
            _df_groups=self._df_groups,
            namespace=self.namespace,
            min_combined_score=int(min_combined_score),
        )

    def mappings(self) -> pl.LazyFrame:
        """Return the input-to-STRING alias mapping lazily.

        Raises:
            CapabilityError: If alias mapping is unavailable for the selected
                namespace or the required aliases source is absent from the
                local snapshot.

        Examples:
            >>> selection.mappings().collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=self._schema_string_map_result,
            batches=lambda request: _iter_mapping_batches(
                snapshot,
                request=request,
            ),
        )

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return normalized input IDs that did not map to STRING lazily.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        if self.namespace == "alias":
            mapped = _mapping_plan(self)[0].select("input_id").unique()
        else:
            if self.dataset.snapshot.file_links is None:
                raise CapabilityError(
                    "STRING links source is absent from this snapshot"
                )
            mapped = pl.concat(
                [
                    self.dataset.scan_links().select(
                        pl.col("protein1").alias("input_id")
                    ),
                    self.dataset.scan_links().select(
                        pl.col("protein2").alias("input_id")
                    ),
                ]
            ).unique()
        input_rows = (
            self._df_group_membership.lazy()
            if self.is_grouped and self._df_group_membership is not None
            else self._df_input_ids.lazy()
        )
        unmatched = input_rows.join(mapped, on="input_id", how="anti")
        return unmatched.select([*self._col_group_id, "input_id"]).sort(
            [*self._col_group_id, "input_id"]
        )

    def edges(self) -> pl.LazyFrame:
        """Return the STRING subnetwork induced by selected IDs lazily.

        Raises:
            CapabilityError: If the required links source is absent from the
                local snapshot, or an alias selection also lacks its required
                aliases source.

        Examples:
            >>> selection.edges().collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=self._schema_edges_result,
            batches=lambda request: _iter_edge_batches(
                snapshot,
                request=request,
            ),
        )


def _iter_mapping_batches(
    selection: StringSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    lf_mapping, schema, _columns = _mapping_plan(selection)
    requested = (
        None
        if request.columns is None
        else [name for name in request.columns if name in schema]
    )
    if requested:
        lf_mapping = lf_mapping.select(requested)
    for frame in lf_mapping.collect_batches(
        chunk_size=request.effective_batch_size,
        engine="streaming",
    ):
        yield frame.cast(
            {name: schema[name] for name in frame.columns if name in schema},
            strict=False,
        )


def _mapping_plan(
    selection: StringSelection,
) -> tuple[pl.LazyFrame, SchemaDict, list[str]]:
    if selection.namespace != "alias":
        raise CapabilityError(
            "STRING alias mapping requires namespace='alias'; "
            "direct namespace='string' selections have no aliases relation"
        )
    alias_info = selection.dataset._alias_schema  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    if alias_info is None:
        raise CapabilityError("STRING aliases source is absent from this snapshot")
    selection.dataset._validate_taxon_compatibility()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    schema = selection._schema_string_map_result  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    if selection._df_input_ids.height == 0:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        return pl.LazyFrame(schema=schema), schema, list(schema)

    lf_aliases = scan_aliases(
        alias_info.file_alias,
        version=selection.dataset.snapshot.parser_version,
    )
    lf_mapping = create_string_mapping_lazy_frame(
        lf_aliases=lf_aliases,
        lf_input_ids=selection._df_input_ids.lazy(),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        source_rank_map=selection.dataset.source_rank_map,
        col_string_id_aliases=alias_info.col_string_id,
        has_source_aliases=alias_info.has_source,
        cols_partition=["input_id", alias_info.col_string_id],
        cols_sort_prefix=["input_id", alias_info.col_string_id],
        cols_select_out=["input_id", alias_info.col_string_id, "alias", "source"],
    )
    if selection.is_grouped and selection._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        columns = [
            "group_id",
            "input_id",
            alias_info.col_string_id,
            "alias",
            "source",
        ]
        lf_mapping = (
            selection._df_group_membership.lazy()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            .join(lf_mapping, on="input_id", how="inner")
            .select(columns)
        )
    else:
        columns = ["input_id", alias_info.col_string_id, "alias", "source"]
        lf_mapping = lf_mapping.select(columns)
    sort_columns = columns
    lf_mapping = lf_mapping.sort(sort_columns)
    return lf_mapping, schema, columns


def _iter_edge_batches(
    selection: StringSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    if selection.dataset.snapshot.file_links is None:
        raise CapabilityError("STRING links source is absent from this snapshot")
    selection.dataset._validate_taxon_compatibility()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    col_group_id = list(selection._col_group_id)  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    if selection.namespace == "string":
        if selection.is_grouped and selection._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            lf_string_ids = selection._df_group_membership.lazy().rename(  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                {"input_id": "string_id"}
            )
        else:
            lf_string_ids = selection._df_input_ids.lazy().rename(  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                {"input_id": "string_id"}
            )
    else:
        lf_mapping, _mapping_schema, _mapping_columns = _mapping_plan(selection)
        alias_info = selection.dataset._alias_schema  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        if alias_info is None:
            raise CapabilityError("STRING aliases source is absent from this snapshot")
        lf_string_ids = (
            lf_mapping.select(  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                [
                    *selection._col_group_id,  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                    alias_info.col_string_id,
                ]
            )
            .rename({alias_info.col_string_id: "string_id"})
            .unique()
        )
    version_link = selection.dataset.snapshot.parser_version
    lf_links = scan_links(selection.dataset.snapshot.file_links, version=version_link)
    validate_required_cols(
        cols_available=lf_links.collect_schema().names(),
        cols_required=SCHEMA_LINKS[version_link].keys(),
        context=f"STRING links file for version {version_link}",
    )
    lf_edges = create_edges_lazy_frame(
        lf_links=lf_links,
        lf_string_ids_a=lf_string_ids.rename({"string_id": "_string_id_a"}),
        lf_string_ids_b=lf_string_ids.rename({"string_id": "_string_id_b"}),
        min_combined_score=selection.min_combined_score,
        cols_join_left_a="_link_a",
        cols_join_right_a="_string_id_a",
        cols_join_left_b=["_link_b", *col_group_id] if col_group_id else "_link_b",
        cols_join_right_b=["_string_id_b", *col_group_id]
        if col_group_id
        else "_string_id_b",
        cols_partition=col_group_id,
        cols_select_out=col_group_id + ["string_id_a", "string_id_b", "combined_score"],
    )
    conflict_keys = [*col_group_id, "_lo", "_hi"]
    conflicts = (
        lf_edges.group_by(conflict_keys)
        .agg(pl.col("combined_score").n_unique().alias("score_count"))
        .filter(pl.col("score_count") > 1)
        .collect(engine="streaming")
    )
    if conflicts.height:
        raise IntegrityError(
            "STRING links contain conflicting combined_score values for "
            f"canonical edges: {conflicts.select(conflict_keys).to_dicts()}"
        )
    lf_output = (
        lf_edges.filter(pl.col("combined_score").ge(int(selection.min_combined_score)))
        .group_by(conflict_keys)
        .agg(pl.col("combined_score").first().cast(pl.Int64))
        .rename({"_lo": "string_id_a", "_hi": "string_id_b"})
        .select(selection._schema_edges_result.keys())  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        .sort(col_group_id + ["string_id_a", "string_id_b"])
    )
    requested = (
        None
        if request.columns is None
        else [
            name
            for name in request.columns
            if name in selection._schema_edges_result  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        ]  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    )
    if requested:
        lf_output = lf_output.select(requested)
    for frame in lf_output.collect_batches(
        chunk_size=request.effective_batch_size,
        engine="streaming",
    ):
        yield frame.cast(
            {
                name: selection._schema_edges_result[name]  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                for name in frame.columns
            },
            strict=False,
        )
