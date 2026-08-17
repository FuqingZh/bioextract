from __future__ import annotations

from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from bioextract._shared import validate_required_cols

from .constant import (
    COLS_MAPPING_RAW,
    COLS_PATHWAY_RAW,
    COLS_RELATION_RAW,
    SCHEMA_MAPPING_RAW,
    SCHEMA_PATHWAY_RAW,
    SCHEMA_RELATION_RAW,
)


def read_mapping_frame(file_uniprot2reactome: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_uniprot2reactome,
        columns=COLS_MAPPING_RAW,
        schema=SCHEMA_MAPPING_RAW,
        context="Reactome UniProt2Reactome file",
    )


def scan_mapping_frame(file_uniprot2reactome: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_uniprot2reactome,
        separator="\t",
        has_header=False,
        new_columns=COLS_MAPPING_RAW,
        schema_overrides=SCHEMA_MAPPING_RAW,
    ).select(COLS_MAPPING_RAW)


def read_pathway_frame(file_pathways: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_pathways,
        columns=COLS_PATHWAY_RAW,
        schema=SCHEMA_PATHWAY_RAW,
        context="Reactome pathways file",
    )


def scan_pathway_frame(file_pathways: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_pathways,
        separator="\t",
        has_header=False,
        new_columns=COLS_PATHWAY_RAW,
        schema_overrides=SCHEMA_PATHWAY_RAW,
    ).select(COLS_PATHWAY_RAW)


def read_relation_frame(file_relations: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_relations,
        columns=COLS_RELATION_RAW,
        schema=SCHEMA_RELATION_RAW,
        context="Reactome pathway relations file",
    )


def scan_relation_frame(file_relations: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_relations,
        separator="\t",
        has_header=False,
        new_columns=COLS_RELATION_RAW,
        schema_overrides=SCHEMA_RELATION_RAW,
    ).select(COLS_RELATION_RAW)


def filter_species_frame(df: pl.DataFrame, species: str | None) -> pl.DataFrame:
    if species is None:
        return df
    return df.filter(pl.col("species") == species)


def filter_relation_frame(
    df_relations: pl.DataFrame,
    df_pathways: pl.DataFrame,
) -> pl.DataFrame:
    df_pathway_ids = df_pathways.select("reactome_pathway_id").unique()
    return (
        df_relations.join(
            df_pathway_ids.rename(
                {"reactome_pathway_id": "parent_reactome_pathway_id"}
            ),
            on="parent_reactome_pathway_id",
            how="inner",
        )
        .join(
            df_pathway_ids.rename({"reactome_pathway_id": "child_reactome_pathway_id"}),
            on="child_reactome_pathway_id",
            how="inner",
        )
        .unique()
        .sort("parent_reactome_pathway_id", "child_reactome_pathway_id")
    )


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    cols_out = [
        "input_id",
        "uniprot_id",
        "reactome_pathway_id",
        "pathway_name",
        "evidence_code",
        "species",
        "reactome_url",
    ]
    df_hits = (
        df_input_ids.join(
            df_mapping,
            left_on="input_id",
            right_on="uniprot_id",
            how="inner",
        )
        .with_columns(pl.col("input_id").alias("uniprot_id"))
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )
    if df_group_membership is None:
        return df_hits
    grouped_cols = ["group_id", *cols_out]
    return (
        df_group_membership.join(df_hits, on="input_id", how="inner")
        .select(grouped_cols)
        .unique()
        .sort(grouped_cols)
    )


def extract_unmatched_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    df_mapped_input_ids = df_mapping.select("input_id").unique().sort("input_id")
    df_unmatched = (
        df_input_ids.join(df_mapped_input_ids, on="input_id", how="anti")
        .select("input_id")
        .sort("input_id")
    )
    if df_group_membership is None:
        return df_unmatched
    return (
        df_group_membership.join(df_unmatched, on="input_id", how="inner")
        .select("group_id", "input_id")
        .sort("group_id", "input_id")
    )


def extract_term2gene_frame(df_mapping: pl.DataFrame) -> pl.DataFrame:
    return (
        df_mapping.select("reactome_pathway_id", "uniprot_id")
        .unique()
        .sort("reactome_pathway_id", "uniprot_id")
    )


def extract_term2name_frame(df_pathways: pl.DataFrame) -> pl.DataFrame:
    return (
        df_pathways.select("reactome_pathway_id", "pathway_name", "species")
        .unique(subset=["reactome_pathway_id"])
        .sort("reactome_pathway_id")
    )


def _read_reactome_tsv(
    file_path: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
    context: str,
) -> pl.DataFrame:
    lf = pl.scan_csv(
        file_path,
        separator="\t",
        has_header=False,
        new_columns=columns,
        schema_overrides=schema,
    )
    df = lf.select(columns).collect()
    validate_required_cols(df.columns, columns, context)
    return df
