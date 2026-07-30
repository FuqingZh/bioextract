from collections.abc import Iterable, Mapping

import polars as pl

from .constant import (
    DEDUP_KEYS_BY_FRAME,
    SCHEMA_ALT_ID,
    SCHEMA_ANCESTOR,
    SCHEMA_DEPTH,
    SCHEMA_EDGE,
    SCHEMA_SUBCELL,
    SCHEMA_SUBSET_DEFINITION,
    SCHEMA_SUBSET_MEMBERSHIP,
    SCHEMA_SYNONYM,
    SCHEMA_TERM,
    SCHEMA_XREF,
)
from .derive import derive_graph_tables
from .model import (
    AltIdColumnBuffer,
    EdgeColumnBuffer,
    SubsetDefinitionColumnBuffer,
    SubsetDefinitionRecord,
    SubsetMembershipColumnBuffer,
    SynonymColumnBuffer,
    TermColumnBuffer,
    TermRecord,
    XrefColumnBuffer,
)
from .parse import (
    normalize_whitespace,
    parse_synonym,
    parse_xref_lossless,
    validate_go_id,
)
from .write import build_deduped_frame


# #region CanonicalColumnBuilders
def generate_alt_id_data_by_cols(go_id: str, alt_ids: list[str]) -> AltIdColumnBuffer:
    alt_id_cols = AltIdColumnBuffer()
    for alt_go_id in alt_ids:
        validate_go_id(alt_go_id)
        if alt_go_id == go_id:
            continue
        alt_id_cols.alt_go_id.append(alt_go_id)
        alt_id_cols.primary_go_id.append(go_id)

    return alt_id_cols


def generate_subset_membership_data_by_cols(
    go_id: str,
    subsets: list[str],
) -> SubsetMembershipColumnBuffer:
    subset_cols = SubsetMembershipColumnBuffer()
    for subset_id in subsets:
        subset_id_normalized = normalize_whitespace(subset_id)
        if not subset_id_normalized:
            continue
        subset_cols.go_id.append(go_id)
        subset_cols.subset_id.append(subset_id_normalized)

    return subset_cols


def generate_subset_definition_data_by_cols(
    records: Iterable[SubsetDefinitionRecord],
) -> SubsetDefinitionColumnBuffer:
    subset_definition_cols = SubsetDefinitionColumnBuffer()
    for record in records:
        subset_definition_cols.subset_id.append(record.subset_id)
        subset_definition_cols.subset_name.append(record.subset_name)

    return subset_definition_cols


def generate_synonym_data_by_cols(
    go_id: str, synonyms: list[str]
) -> SynonymColumnBuffer:
    synonym_cols = SynonymColumnBuffer()
    for raw_synonym in synonyms:
        synonym_record = parse_synonym(raw_synonym)
        if synonym_record is None:
            continue

        synonym_cols.go_id.append(go_id)
        synonym_cols.synonym_text.append(synonym_record.synonym_text)
        synonym_cols.synonym_scope.append(synonym_record.synonym_scope)
        synonym_cols.synonym_type_name.append(synonym_record.synonym_type_name)
        synonym_cols.dbxref_text.append(synonym_record.dbxref_text)

    return synonym_cols


def generate_xref_data_by_cols(go_id: str, xrefs: list[str]) -> XrefColumnBuffer:
    xref_cols = XrefColumnBuffer()
    for raw_xref in xrefs:
        xref_text = normalize_whitespace(raw_xref)
        xref_parts = parse_xref_lossless(xref_text)
        xref_cols.go_id.append(go_id)
        xref_cols.xref_text.append(xref_text)
        xref_cols.xref_db.append(xref_parts.xref_db)
        xref_cols.xref_id.append(xref_parts.xref_id)

    return xref_cols


# #endregion
################################################################################
# #region FrameBuilders
def _build_tidy_frames(
    records: Iterable[TermRecord],
    *,
    subset_definitions: Iterable[SubsetDefinitionRecord] = (),
) -> dict[str, pl.DataFrame]:
    """Build canonical and graph-derived GO ontology frames.

    Args:
        records: Parsed GO term records. The iterable is consumed once.
        subset_definitions: Optional OBO header definitions used to populate
            the subset lookup frame.

    Returns:
        DataFrames keyed by `term`, `edge`, `synonym`, `xref`, `alt_id`,
        `subset_membership`, `subset_definition`, `ancestor_all`, and `depth`.
        Every frame is deduplicated according to its canonical key.

    Raises:
        ValueError: If an alternate GO identifier is invalid or hierarchical
            edges contain a cycle, which prevents ancestor and depth
            derivation.

    Notes:
        This builder consumes parser-owned record types. Callers should use
        `GODatabase.build_tidy()` rather than depend on those internal models.
    """
    term_data = TermColumnBuffer()
    edge_data = EdgeColumnBuffer()
    alt_id_data = AltIdColumnBuffer()
    subset_membership_data = SubsetMembershipColumnBuffer()
    synonym_data = SynonymColumnBuffer()
    xref_data = XrefColumnBuffer()

    for record in records:
        term_data.append_record(record)
        for edge in record.parents:
            edge_data.append_edge(edge)
        alt_id_data.extend(generate_alt_id_data_by_cols(record.go_id, record.alt_ids))
        subset_membership_data.extend(
            generate_subset_membership_data_by_cols(record.go_id, record.subsets)
        )
        synonym_data.extend(
            generate_synonym_data_by_cols(record.go_id, record.synonyms)
        )
        xref_data.extend(generate_xref_data_by_cols(record.go_id, record.xrefs))

    df_term = build_deduped_frame(
        term_data,
        schema=SCHEMA_TERM,
        dedup_keys=DEDUP_KEYS_BY_FRAME["term"],
    )
    df_edge = build_deduped_frame(
        edge_data,
        schema=SCHEMA_EDGE,
        dedup_keys=DEDUP_KEYS_BY_FRAME["edge"],
    )
    df_synonym = build_deduped_frame(
        synonym_data,
        schema=SCHEMA_SYNONYM,
        dedup_keys=DEDUP_KEYS_BY_FRAME["synonym"],
    )
    df_xref = build_deduped_frame(
        xref_data,
        schema=SCHEMA_XREF,
        dedup_keys=DEDUP_KEYS_BY_FRAME["xref"],
    )
    df_alt_id = build_deduped_frame(
        alt_id_data,
        schema=SCHEMA_ALT_ID,
        dedup_keys=DEDUP_KEYS_BY_FRAME["alt_id"],
    )
    df_subset_membership = build_deduped_frame(
        subset_membership_data,
        schema=SCHEMA_SUBSET_MEMBERSHIP,
        dedup_keys=DEDUP_KEYS_BY_FRAME["subset_membership"],
    )
    df_subset_definition = build_deduped_frame(
        generate_subset_definition_data_by_cols(subset_definitions),
        schema=SCHEMA_SUBSET_DEFINITION,
        dedup_keys=DEDUP_KEYS_BY_FRAME["subset_definition"],
    )

    ancestor_data, depth_data = derive_graph_tables(term_data, edge_data)
    df_ancestor = build_deduped_frame(
        ancestor_data,
        schema=SCHEMA_ANCESTOR,
        dedup_keys=("go_id", "ancestor_go_id"),
    )
    df_depth = build_deduped_frame(
        depth_data,
        schema=SCHEMA_DEPTH,
        dedup_keys=("go_id",),
    )

    return {
        "term": df_term,
        "edge": df_edge,
        "synonym": df_synonym,
        "xref": df_xref,
        "alt_id": df_alt_id,
        "subset_membership": df_subset_membership,
        "subset_definition": df_subset_definition,
        "ancestor_all": df_ancestor,
        "depth": df_depth,
    }


def extract_subcell_frame(
    frames: Mapping[str, pl.DataFrame],
    *,
    include_obsolete: bool = False,
) -> pl.DataFrame:
    """Project GO cellular-component terms into the subcell output schema.

    Args:
        frames: Tidy frame mapping containing `term` and `depth`.
        include_obsolete: Whether to retain obsolete cellular-component terms.

    Returns:
        A table sorted by `go_id` with subcell names, definitions, and minimum
        and maximum depths from an ontology root.

    Raises:
        KeyError: If either required frame is absent.
    """
    df_term = frames["term"]
    df_depth = frames["depth"]
    df_subcell = (
        df_term.filter(pl.col("namespace") == "cellular_component")
        .join(
            df_depth.select("go_id", "min_depth_from_root", "max_depth_from_root"),
            on="go_id",
            how="left",
        )
        .select(
            "go_id",
            pl.col("term_name").alias("subcell_name"),
            "definition",
            "is_obsolete",
            "min_depth_from_root",
            "max_depth_from_root",
        )
    )
    if not include_obsolete:
        df_subcell = df_subcell.filter(~pl.col("is_obsolete"))
    return df_subcell.select(SCHEMA_SUBCELL.keys()).sort("go_id")


# #endregion
