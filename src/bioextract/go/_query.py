from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import register_replayable_source

if TYPE_CHECKING:
    from .go import GODatabase


GO_UNMATCHED_REASONS = (
    "not_found",
    "obsolete_excluded",
    "no_matching_ancestor",
)
GO_UNMATCHED_REASON_DTYPE = pl.Enum(GO_UNMATCHED_REASONS)

GO_ANCESTOR_MATCH_SCHEMA: SchemaDict = {
    "input_go_id": pl.String,
    "go_id": pl.String,
    "ancestor_go_id": pl.String,
    "ancestor_term_name": pl.String,
    "ancestor_namespace": pl.String,
    "min_distance": pl.Int32,
    "target_subset_id": pl.String,
}
GO_ANCESTOR_UNMATCHED_SCHEMA: SchemaDict = {
    "input_go_id": pl.String,
    "resolved_go_id": pl.String,
    "reason": GO_UNMATCHED_REASON_DTYPE,
}

_GO_QUERY_RESULT_SCHEMA: SchemaDict = {
    "result_kind": pl.String,
    "input_order": pl.Int64,
    "input_go_id": pl.String,
    "go_id": pl.String,
    "ancestor_go_id": pl.String,
    "ancestor_term_name": pl.String,
    "ancestor_namespace": pl.String,
    "min_distance": pl.Int32,
    "reason": pl.String,
}

_GO_TERM_SCHEMA: SchemaDict = {
    "go_id": pl.String,
    "term_name": pl.String,
    "namespace": pl.String,
    "definition": pl.String,
    "is_obsolete": pl.Boolean,
    "comment": pl.String,
}


def empty_ancestor_matches() -> pl.DataFrame:
    return pl.DataFrame(schema=GO_ANCESTOR_MATCH_SCHEMA)


def empty_ancestor_unmatched() -> pl.DataFrame:
    return pl.DataFrame(schema=GO_ANCESTOR_UNMATCHED_SCHEMA)


def select_terms_from_publication(
    database: GODatabase,
    *,
    df_input_terms: pl.DataFrame | None,
    namespace: str | None,
    subset_id: str | None,
    include_obsolete: bool,
    resolve_alt_ids: bool,
) -> pl.DataFrame:
    """Select GO terms with publication-side filtering and joins."""
    if df_input_terms is not None and df_input_terms.is_empty():
        database._assert_publication_identity()  # pyright: ignore[reportPrivateUsage]
        schema: SchemaDict = dict(_GO_TERM_SCHEMA)
        schema = {
            **{"input_go_id": pl.String},
            **schema,
            **({"subset_id": pl.String} if subset_id is not None else {}),
        }
        return pl.DataFrame(schema=schema)

    params: list[object] = []
    ctes: list[str] = []
    order_by = "term.go_id"
    select_prefix = ""
    join_clause = ""

    if df_input_terms is not None:
        input_ids = df_input_terms.get_column("input_go_id").to_list()
        params.append(input_ids)
        ctes.extend(
            [
                "input AS ("
                " SELECT input_go_id, input_order - 1 AS input_order"
                " FROM unnest(?) WITH ORDINALITY AS input(input_go_id, input_order)"
                ")",
                "resolved_ids AS ("
                " SELECT input.input_go_id, input.input_order,"
                " CASE WHEN primary_term.go_id IS NOT NULL"
                " THEN primary_term.go_id"
                + (
                    " ELSE alternate.primary_go_id"
                    if resolve_alt_ids
                    else " ELSE NULL::VARCHAR"
                )
                + " END AS go_id"
                " FROM input"
                " LEFT JOIN term AS primary_term"
                " ON primary_term.go_id = input.input_go_id"
                + (
                    " LEFT JOIN term_alternate_id AS alternate"
                    " ON alternate.alt_go_id = input.input_go_id"
                    if resolve_alt_ids
                    else ""
                )
                + ")",
            ]
        )
        select_prefix = "resolved.input_go_id, "
        join_clause = "JOIN resolved_ids AS resolved ON resolved.go_id = term.go_id"
        order_by = "resolved.input_order, term.go_id"

    where_clauses: list[str] = []
    if not include_obsolete:
        where_clauses.append("NOT term.is_obsolete")
    if namespace is not None:
        where_clauses.append("term.namespace = ?")
        params.append(namespace)
    if subset_id is not None:
        where_clauses.append(
            "EXISTS (SELECT 1 FROM subset_membership AS membership "
            "WHERE membership.go_id = term.go_id AND membership.subset_id = ?)"
        )
        params.append(subset_id)

    where_sql = "" if not where_clauses else "WHERE " + " AND ".join(where_clauses)
    prefix_sql = "" if not ctes else "WITH " + ", ".join(ctes) + "\n"

    query = (
        prefix_sql
        + "SELECT "
        + select_prefix
        + "term.go_id, term.term_name, term.namespace, term.definition, "
        + "term.is_obsolete, term.comment "
        + "FROM term "
        + join_clause
        + " "
        + where_sql
        + " ORDER BY "
        + order_by
    )

    with database.connect() as connection:
        rows = connection.execute(query, params).fetchall()

    schema: SchemaDict = dict(_GO_TERM_SCHEMA)
    if df_input_terms is not None:
        schema = {"input_go_id": pl.String, **schema}
    if not rows:
        frame = pl.DataFrame(schema=schema)
    else:
        frame = pl.DataFrame(rows, schema=schema, orient="row")
    if subset_id is not None:
        frame = frame.with_columns(pl.lit(subset_id).cast(pl.String).alias("subset_id"))
    return frame


def list_subsets_from_publication(database: GODatabase) -> pl.DataFrame:
    """List GO subsets with counts computed inside publication DuckDB."""
    query = """
        WITH counts AS (
            SELECT subset_id, count(DISTINCT go_id)::UBIGINT AS num_terms
            FROM subset_membership
            GROUP BY subset_id
        )
        SELECT
            COALESCE(definition.subset_id, counts.subset_id) AS subset_id,
            definition.subset_name,
            COALESCE(counts.num_terms, 0)::UBIGINT AS num_terms
        FROM subset_definition AS definition
        FULL OUTER JOIN counts USING (subset_id)
        ORDER BY subset_id
    """
    with database.connect() as connection:
        rows = connection.execute(query).fetchall()
    schema: SchemaDict = {
        "subset_id": pl.String,
        "subset_name": pl.String,
        "num_terms": pl.UInt32,
    }
    if not rows:
        return pl.DataFrame(schema=schema)
    return pl.DataFrame(rows, schema=schema, orient="row")


@dataclass(slots=True)
class GOAncestorSelection:
    """Deferred selection of GO ontology ancestors.

    The selection supports both an OBO source handle and a validated DuckDB
    publication. Its terminals retain input lineage and share one selection-
    local resolution cache.

    Examples:
        Extract projected ancestor rows from a selection:

        >>> selection.ancestors().collect_schema().names()  # doctest: +SKIP
        ['input_go_id', 'go_id', 'ancestor_go_id', 'ancestor_term_name', 'ancestor_namespace', 'min_distance', 'target_subset_id']
    """

    database: GODatabase = field(repr=False)
    _df_input_terms: pl.DataFrame = field(repr=False)
    target_subset_id: str | None
    include_self: bool
    resolve_alt_ids: bool
    include_obsolete: bool
    _df_matches: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_unmatched: pl.DataFrame | None = field(default=None, init=False, repr=False)

    def ancestors(self) -> pl.LazyFrame:
        """Return accepted ancestor and optional self mappings lazily.

        Examples:
            >>> selection.ancestors().collect()  # doctest: +SKIP
            shape: (..., 7)
        """
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=GO_ANCESTOR_MATCH_SCHEMA,
            batches=lambda request: _iter_single_frame(
                snapshot._eager_ancestors(), request.effective_batch_size
            ),
        )

    def _eager_ancestors(self) -> pl.DataFrame:
        """Return accepted ancestor and optional self mappings.

        Examples:
            Read the canonical ancestor ID from a selection:

            >>> selection.ancestors().collect().get_column("ancestor_go_id")  # doctest: +SKIP
            shape: (1,)
            Series: 'ancestor_go_id' [str]
            ["GO:0000001"]
        """
        self._resolve()
        assert self._df_matches is not None
        return self._df_matches

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return inputs without an accepted ancestor lazily.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=GO_ANCESTOR_UNMATCHED_SCHEMA,
            batches=lambda request: _iter_single_frame(
                snapshot._eager_unmatched_ids(), request.effective_batch_size
            ),
        )

    def _eager_unmatched_ids(self) -> pl.DataFrame:
        """Return inputs without an accepted ancestor and their reason.

        Examples:
            Inspect the reason for an absent input:

            >>> selection.unmatched_ids().collect_schema().names()  # doctest: +SKIP
            ['input_go_id', 'resolved_go_id', 'reason']
        """
        self._resolve()
        assert self._df_unmatched is not None
        return self._df_unmatched

    def _resolve(self) -> None:
        if self._df_matches is not None and self._df_unmatched is not None:
            return
        if self._df_input_terms.is_empty():
            if self.database._publication_path is not None:  # pyright: ignore[reportPrivateUsage]
                self.database._assert_publication_identity()  # pyright: ignore[reportPrivateUsage]
            self._df_matches = empty_ancestor_matches()
            self._df_unmatched = empty_ancestor_unmatched()
            return
        if self.database._publication_path is not None:  # pyright: ignore[reportPrivateUsage]
            df_matches, df_unmatched = _resolve_publication(self)
        else:
            df_matches, df_unmatched = _resolve_source(self)
        self._df_matches = df_matches
        self._df_unmatched = df_unmatched


def _iter_single_frame(frame: pl.DataFrame, batch_size: int):
    for offset in range(0, frame.height, batch_size):
        yield frame.slice(offset, batch_size)


def _resolve_publication(
    selection: GOAncestorSelection,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    input_ids = selection._df_input_terms.get_column(  # pyright: ignore[reportPrivateUsage]
        "input_go_id"
    ).to_list()
    params: list[object] = [input_ids]
    alt_join = ""
    alt_expr = "NULL::VARCHAR"
    if selection.resolve_alt_ids:
        alt_join = (
            " LEFT JOIN term_alternate_id AS alternate"
            " ON alternate.alt_go_id = input.input_go_id"
        )
        alt_expr = "alternate.primary_go_id"

    input_policy = "TRUE" if selection.include_obsolete else "NOT resolved.is_obsolete"
    ancestor_policy = (
        "TRUE" if selection.include_obsolete else "NOT ancestor_term.is_obsolete"
    )
    self_sql = ""
    if selection.include_self:
        self_sql = (
            " UNION ALL "
            "SELECT resolved.input_order, resolved.input_go_id, resolved.go_id,"
            " resolved.go_id AS ancestor_go_id, CAST(0 AS INTEGER) AS min_distance"
            " FROM resolved"
            " WHERE resolved.go_id IS NOT NULL AND " + input_policy
        )

    subset_sql = ""
    if selection.target_subset_id is not None:
        subset_sql = (
            " AND EXISTS (SELECT 1 FROM subset_membership AS membership"
            " WHERE membership.go_id = candidate.ancestor_go_id"
            " AND membership.subset_id = ?)"
        )
        params.append(selection.target_subset_id)

    obsolete_reason_sql = (
        " WHEN resolved.is_obsolete THEN 'obsolete_excluded'"
        if not selection.include_obsolete
        else ""
    )

    query = (
        "WITH input AS ("
        " SELECT input_go_id, input_order - 1 AS input_order"
        " FROM unnest(?) WITH ORDINALITY AS input(input_go_id, input_order)"
        "), resolved_ids AS ("
        " SELECT input.input_order, input.input_go_id,"
        " CASE WHEN primary_term.go_id IS NOT NULL THEN primary_term.go_id"
        f" ELSE {alt_expr} END AS go_id"
        " FROM input"
        " LEFT JOIN term AS primary_term"
        " ON primary_term.go_id = input.input_go_id"
        f"{alt_join}"
        "), resolved AS ("
        " SELECT resolved_ids.*, term.is_obsolete"
        " FROM resolved_ids"
        " LEFT JOIN term ON term.go_id = resolved_ids.go_id"
        "), candidate_rows AS ("
        " SELECT resolved.input_order, resolved.input_go_id, resolved.go_id,"
        " ancestor.ancestor_go_id, ancestor.min_distance"
        " FROM resolved"
        " JOIN term_ancestor AS ancestor ON ancestor.go_id = resolved.go_id"
        " JOIN term AS ancestor_term ON ancestor_term.go_id = ancestor.ancestor_go_id"
        " WHERE resolved.go_id IS NOT NULL AND "
        + input_policy
        + self_sql
        + "), accepted_rows AS ("
        " SELECT candidate.input_order, candidate.input_go_id, candidate.go_id,"
        " candidate.ancestor_go_id, candidate.min_distance,"
        " ancestor_term.term_name AS ancestor_term_name,"
        " ancestor_term.namespace AS ancestor_namespace"
        " FROM candidate_rows AS candidate"
        " JOIN term AS ancestor_term"
        " ON ancestor_term.go_id = candidate.ancestor_go_id"
        " WHERE "
        + ancestor_policy.replace("ancestor_term", "ancestor_term")
        + subset_sql
        + "), ranked_matches AS ("
        " SELECT accepted_rows.*, ROW_NUMBER() OVER ("
        " PARTITION BY input_go_id, go_id, ancestor_go_id"
        " ORDER BY min_distance"
        " ) AS row_number"
        " FROM accepted_rows"
        "), matches AS ("
        " SELECT input_order, input_go_id, go_id, ancestor_go_id, min_distance,"
        " ancestor_term_name, ancestor_namespace"
        " FROM ranked_matches WHERE row_number = 1"
        "), matched_inputs AS ("
        " SELECT DISTINCT input_go_id FROM matches"
        "), combined AS ("
        " SELECT 'match' AS result_kind, input_order, input_go_id, go_id,"
        " ancestor_go_id, ancestor_term_name, ancestor_namespace, min_distance,"
        " NULL::VARCHAR AS reason FROM matches"
        " UNION ALL "
        " SELECT 'unmatched' AS result_kind, resolved.input_order,"
        " resolved.input_go_id, resolved.go_id, NULL::VARCHAR, NULL::VARCHAR,"
        " NULL::VARCHAR, NULL::INTEGER,"
        " CASE WHEN resolved.go_id IS NULL THEN 'not_found'"
        + obsolete_reason_sql
        + " ELSE 'no_matching_ancestor' END AS reason"
        " FROM resolved"
        " LEFT JOIN matched_inputs USING (input_go_id)"
        " WHERE matched_inputs.input_go_id IS NULL"
        ") SELECT * FROM combined"
        " ORDER BY input_order, result_kind, min_distance NULLS LAST, ancestor_go_id NULLS LAST"
    )

    with selection.database.connect() as connection:
        rows = connection.execute(query, params).fetchall()
    combined = pl.DataFrame(rows, schema=_GO_QUERY_RESULT_SCHEMA, orient="row")
    return _split_query_results(combined, target_subset_id=selection.target_subset_id)


def _resolve_source(
    selection: GOAncestorSelection,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    frames = selection.database._collect_frames(  # pyright: ignore[reportPrivateUsage]
        {"term", "alt_id", "ancestor_all", "subset_membership"}
    )
    df_input = selection._df_input_terms  # pyright: ignore[reportPrivateUsage]
    df_term = frames["term"]
    df_resolved = df_input.join(
        df_term.select(
            pl.col("go_id").alias("direct_lookup_go_id"),
            pl.col("go_id").alias("direct_go_id"),
        ),
        left_on="input_go_id",
        right_on="direct_lookup_go_id",
        how="left",
    )
    if selection.resolve_alt_ids:
        df_resolved = df_resolved.join(
            frames["alt_id"],
            left_on="input_go_id",
            right_on="alt_go_id",
            how="left",
        )
        df_resolved = df_resolved.with_columns(
            pl.coalesce("direct_go_id", "primary_go_id").alias("go_id")
        )
    else:
        df_resolved = df_resolved.with_columns(pl.col("direct_go_id").alias("go_id"))
    df_resolved = df_resolved.join(
        df_term.select("go_id", "is_obsolete"), on="go_id", how="left"
    ).select("input_go_id", "input_order", "go_id", "is_obsolete")

    df_candidates = (
        df_resolved.filter(
            pl.col("go_id").is_not_null()
            & (pl.lit(selection.include_obsolete) | ~pl.col("is_obsolete"))
        )
        .join(frames["ancestor_all"], on="go_id", how="inner")
        .select(
            "input_go_id",
            "input_order",
            "go_id",
            pl.col("ancestor_go_id"),
            "min_distance",
        )
    )
    if selection.include_self:
        df_self = df_resolved.filter(
            pl.col("go_id").is_not_null()
            & (pl.lit(selection.include_obsolete) | ~pl.col("is_obsolete"))
        ).select(
            "input_go_id",
            "input_order",
            "go_id",
            pl.col("go_id").alias("ancestor_go_id"),
            pl.lit(0, dtype=pl.Int32).alias("min_distance"),
        )
        df_candidates = pl.concat([df_candidates, df_self], how="vertical")

    df_candidates = df_candidates.join(
        df_term.select(
            pl.col("go_id").alias("ancestor_go_id"),
            pl.col("term_name").alias("ancestor_term_name"),
            pl.col("namespace").alias("ancestor_namespace"),
            pl.col("is_obsolete").alias("ancestor_is_obsolete"),
        ),
        on="ancestor_go_id",
        how="inner",
    )
    if not selection.include_obsolete:
        df_candidates = df_candidates.filter(~pl.col("ancestor_is_obsolete"))
    if selection.target_subset_id is not None:
        df_subset = frames["subset_membership"].filter(
            pl.col("subset_id") == selection.target_subset_id
        )
        df_candidates = df_candidates.join(
            df_subset.select(pl.col("go_id").alias("ancestor_go_id")),
            on="ancestor_go_id",
            how="inner",
        )

    df_matches = (
        df_candidates.sort("input_order", "min_distance", "ancestor_go_id")
        .unique(
            subset=["input_go_id", "go_id", "ancestor_go_id"],
            keep="first",
            maintain_order=True,
        )
        .with_columns(
            pl.lit(selection.target_subset_id, dtype=pl.String).alias(
                "target_subset_id"
            )
        )
        .select(
            "input_go_id",
            "go_id",
            "ancestor_go_id",
            "ancestor_term_name",
            "ancestor_namespace",
            "min_distance",
            "target_subset_id",
        )
    )
    matched_inputs = df_matches.select("input_go_id").unique()
    df_unmatched = (
        df_resolved.join(matched_inputs, on="input_go_id", how="anti")
        .with_columns(
            pl.when(pl.col("go_id").is_null())
            .then(pl.lit("not_found"))
            .when(~pl.lit(selection.include_obsolete) & pl.col("is_obsolete"))
            .then(pl.lit("obsolete_excluded"))
            .otherwise(pl.lit("no_matching_ancestor"))
            .cast(GO_UNMATCHED_REASON_DTYPE)
            .alias("reason")
        )
        .select(
            "input_go_id",
            pl.col("go_id").alias("resolved_go_id"),
            "reason",
        )
    )
    return df_matches, df_unmatched


def _split_query_results(
    combined: pl.DataFrame,
    *,
    target_subset_id: str | None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    if combined.is_empty():
        return empty_ancestor_matches(), empty_ancestor_unmatched()

    df_matches = (
        combined.filter(pl.col("result_kind") == "match")
        .with_columns(
            pl.lit(target_subset_id, dtype=pl.String).alias("target_subset_id")
        )
        .select(
            "input_go_id",
            "go_id",
            "ancestor_go_id",
            "ancestor_term_name",
            "ancestor_namespace",
            "min_distance",
            "target_subset_id",
        )
    )
    df_unmatched = combined.filter(pl.col("result_kind") == "unmatched").select(
        "input_go_id",
        pl.col("go_id").alias("resolved_go_id"),
        pl.col("reason").cast(GO_UNMATCHED_REASON_DTYPE).alias("reason"),
    )
    return df_matches, df_unmatched
