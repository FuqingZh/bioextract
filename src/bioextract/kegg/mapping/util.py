from __future__ import annotations

from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from .constant import (
    COLS_MAPPING,
    KIND_INPUT_ID_VALUES,
    SCHEMA_MAPPING,
    KeggInputIdKind,
)


def read_conv_uniprot_frame(file_conv_uniprot: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_conv_uniprot,
        columns=["UniProtId", "KeggGeneId"],
        schema={"UniProtId": pl.String, "KeggGeneId": pl.String},
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("UniProtId").str.strip_prefix("up:"))
        .unique()
        .sort("KeggGeneId", "UniProtId")
    )


def read_conv_ncbi_geneid_frame(file_conv_ncbi_geneid: Path | None) -> pl.DataFrame:
    schema = {"NcbiGeneId": pl.String, "KeggGeneId": pl.String}
    if file_conv_ncbi_geneid is None:
        return pl.DataFrame(schema=schema)
    df = _read_two_col_tsv(
        file_conv_ncbi_geneid,
        columns=["NcbiGeneId", "KeggGeneId"],
        schema=schema,
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("NcbiGeneId").str.strip_prefix("ncbi-geneid:"))
        .unique()
        .sort("KeggGeneId", "NcbiGeneId")
    )


def read_gene_ko_frame(file_gene_ko: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_gene_ko,
        columns=["KeggGeneId", "KoId"],
        schema={"KeggGeneId": pl.String, "KoId": pl.String},
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("KoId").str.strip_prefix("ko:"))
        .unique()
        .sort("KeggGeneId", "KoId")
    )


def read_gene_pathway_frame(file_gene_pathway: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_gene_pathway,
        columns=["KeggGeneId", "KeggPathwayId"],
        schema={"KeggGeneId": pl.String, "KeggPathwayId": pl.String},
    )
    if df.height == 0:
        return df.with_columns(pl.lit(None, dtype=pl.String).alias("PathwayMapId"))
    return (
        df.with_columns(pl.col("KeggPathwayId").str.strip_prefix("path:"))
        .with_columns(derive_pathway_map_id_expr().alias("PathwayMapId"))
        .unique()
        .sort("KeggGeneId", "KeggPathwayId")
    )


def read_gene_list_frame(file_gene_list: Path | None) -> pl.DataFrame:
    schema = {
        "KeggGeneId": pl.String,
        "GeneSymbol": pl.String,
        "GeneDescription": pl.String,
    }
    if file_gene_list is None:
        return pl.DataFrame(schema=schema)
    df = _read_two_col_tsv(
        file_gene_list,
        columns=["KeggGeneId", "_GeneText"],
        schema={"KeggGeneId": pl.String, "_GeneText": pl.String},
    )
    if df.height == 0:
        return pl.DataFrame(schema=schema)

    rows: list[dict[str, str | None]] = []
    for row in df.iter_rows(named=True):
        text = str(row["_GeneText"]).strip()
        if ";" in text:
            symbol, description = text.split(";", 1)
            symbol = symbol.strip() or None
            description = description.strip() or None
        else:
            symbol = text or None
            description = None
        rows.append(
            {
                "KeggGeneId": row["KeggGeneId"],
                "GeneSymbol": symbol,
                "GeneDescription": description,
            }
        )
    return pl.DataFrame(rows, schema=schema).unique().sort("KeggGeneId")


def build_mapping_frame(
    *,
    organism_code: str,
    df_conv_uniprot: pl.DataFrame,
    df_conv_ncbi_geneid: pl.DataFrame,
    df_gene_ko: pl.DataFrame,
    df_gene_pathway: pl.DataFrame,
    df_gene_list: pl.DataFrame,
) -> pl.DataFrame:
    df_gene_ids = (
        pl.concat(
            [
                df_conv_uniprot.select("KeggGeneId"),
                df_conv_ncbi_geneid.select("KeggGeneId"),
                df_gene_ko.select("KeggGeneId"),
                df_gene_pathway.select("KeggGeneId"),
                df_gene_list.select("KeggGeneId"),
            ],
            how="vertical_relaxed",
        )
        .drop_nulls("KeggGeneId")
        .unique()
    )
    if df_gene_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_MAPPING)

    validate_organism_code(df_gene_ids, organism_code=organism_code)
    return (
        df_gene_ids.join(df_conv_uniprot, on="KeggGeneId", how="left")
        .join(df_conv_ncbi_geneid, on="KeggGeneId", how="left")
        .join(df_gene_ko, on="KeggGeneId", how="left")
        .join(df_gene_pathway, on="KeggGeneId", how="left")
        .join(df_gene_list, on="KeggGeneId", how="left")
        .with_columns(pl.lit(organism_code).alias("OrganismCode"))
        .select(COLS_MAPPING)
        .unique()
        .sort(COLS_MAPPING)
    )


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    kind_input_id: KeggInputIdKind,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    validate_kind_input_id(kind_input_id)
    col_join = col_join_by_kind(kind_input_id)
    cols_group = list(cols_group_id)
    cols_out = cols_group + ["InputId", "KindInputId"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="InputId",
            right_on=col_join,
            how="inner",
        )
        .with_columns(
            pl.col("InputId").alias(col_join),
            pl.lit(kind_input_id).alias("KindInputId"),
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


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


def validate_kind_input_id(kind_input_id: str) -> None:
    if kind_input_id not in KIND_INPUT_ID_VALUES:
        raise ValueError(
            "kind_input_id must be one of: "
            f"{', '.join(KIND_INPUT_ID_VALUES)}; got {kind_input_id!r}"
        )


def col_join_by_kind(kind_input_id: KeggInputIdKind) -> str:
    if kind_input_id == "uniprot":
        return "UniProtId"
    if kind_input_id == "ncbi_geneid":
        return "NcbiGeneId"
    if kind_input_id == "kegg_gene":
        return "KeggGeneId"
    validate_kind_input_id(kind_input_id)
    raise AssertionError("unreachable")


def validate_organism_code(df_gene_ids: pl.DataFrame, *, organism_code: str) -> None:
    prefix = f"{organism_code}:"
    df_invalid = df_gene_ids.filter(~pl.col("KeggGeneId").str.starts_with(prefix))
    if df_invalid.height == 0:
        return
    examples = df_invalid.get_column("KeggGeneId").head(5).to_list()
    raise ValueError(
        "KEGG gene IDs do not match organism_code: "
        f"organism_code={organism_code!r}, examples={examples!r}"
    )


def derive_pathway_map_id_expr() -> pl.Expr:
    return pl.lit("map") + pl.col("KeggPathwayId").str.extract(
        r"^[A-Za-z]+([0-9]{5})$", 1
    )


def _read_two_col_tsv(
    file_path: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
) -> pl.DataFrame:
    if file_path.stat().st_size == 0:
        return pl.DataFrame(schema=schema)
    return (
        pl.scan_csv(
            file_path,
            separator="\t",
            has_header=False,
            new_columns=columns,
            schema_overrides=schema,
        )
        .select(columns)
        .collect()
    )
