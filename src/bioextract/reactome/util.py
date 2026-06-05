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


def read_pathway_frame(file_pathways: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_pathways,
        columns=COLS_PATHWAY_RAW,
        schema=SCHEMA_PATHWAY_RAW,
        context="Reactome pathways file",
    )


def read_relation_frame(file_relations: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_relations,
        columns=COLS_RELATION_RAW,
        schema=SCHEMA_RELATION_RAW,
        context="Reactome pathway relations file",
    )


def filter_species_frame(df: pl.DataFrame, species: str | None) -> pl.DataFrame:
    if species is None:
        return df
    return df.filter(pl.col("Species") == species)


def filter_relation_frame(
    df_relations: pl.DataFrame,
    df_pathways: pl.DataFrame,
) -> pl.DataFrame:
    df_pathway_ids = df_pathways.select("ReactomePathwayId").unique()
    return (
        df_relations.join(
            df_pathway_ids.rename({"ReactomePathwayId": "ParentReactomePathwayId"}),
            on="ParentReactomePathwayId",
            how="inner",
        )
        .join(
            df_pathway_ids.rename({"ReactomePathwayId": "ChildReactomePathwayId"}),
            on="ChildReactomePathwayId",
            how="inner",
        )
        .unique()
        .sort("ParentReactomePathwayId", "ChildReactomePathwayId")
    )


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    cols_group = list(cols_group_id)
    cols_out = cols_group + [
        "InputId",
        "UniProtId",
        "ReactomePathwayId",
        "PathwayName",
        "EvidenceCode",
        "Species",
        "ReactomeUrl",
    ]
    df_hits = (
        df_input_ids.join(
            df_mapping,
            left_on="InputId",
            right_on="UniProtId",
            how="inner",
        )
        .with_columns(pl.col("InputId").alias("UniProtId"))
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )
    return df_hits


def extract_unmapped_input_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    cols_index = list(cols_group_id) + ["InputId"]
    df_mapped_input_ids = df_mapping.select(cols_index).unique().sort(cols_index)
    return (
        df_input_ids.join(df_mapped_input_ids, on=cols_index, how="anti")
        .select(cols_index)
        .sort(cols_index)
    )


def extract_term2gene_frame(df_mapping: pl.DataFrame) -> pl.DataFrame:
    return (
        df_mapping.select("ReactomePathwayId", "UniProtId")
        .unique()
        .sort("ReactomePathwayId", "UniProtId")
    )


def extract_term2name_frame(df_pathways: pl.DataFrame) -> pl.DataFrame:
    return (
        df_pathways.select("ReactomePathwayId", "PathwayName", "Species")
        .unique(subset=["ReactomePathwayId"])
        .sort("ReactomePathwayId")
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
