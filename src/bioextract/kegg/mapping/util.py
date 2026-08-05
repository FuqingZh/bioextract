from __future__ import annotations

from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from .constant import (
    COLS_MAPPING,
    NAMESPACE_VALUES,
    SCHEMA_MAPPING,
    KEGGNamespace,
)


def read_conv_uniprot_frame(file_conv_uniprot: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_conv_uniprot,
        columns=["uniprot_id", "kegg_gene_id"],
        schema={"uniprot_id": pl.String, "kegg_gene_id": pl.String},
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("uniprot_id").str.strip_prefix("up:"))
        .unique()
        .sort("kegg_gene_id", "uniprot_id")
    )


def read_conv_ncbi_geneid_frame(file_conv_ncbi_geneid: Path | None) -> pl.DataFrame:
    schema = {"ncbi_gene_id": pl.String, "kegg_gene_id": pl.String}
    if file_conv_ncbi_geneid is None:
        return pl.DataFrame(schema=schema)
    df = _read_two_col_tsv(
        file_conv_ncbi_geneid,
        columns=["ncbi_gene_id", "kegg_gene_id"],
        schema=schema,
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("ncbi_gene_id").str.strip_prefix("ncbi-geneid:"))
        .unique()
        .sort("kegg_gene_id", "ncbi_gene_id")
    )


def read_gene_ko_frame(file_gene_ko: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_gene_ko,
        columns=["kegg_gene_id", "ko_id"],
        schema={"kegg_gene_id": pl.String, "ko_id": pl.String},
    )
    if df.height == 0:
        return df
    return (
        df.with_columns(pl.col("ko_id").str.strip_prefix("ko:"))
        .unique()
        .sort("kegg_gene_id", "ko_id")
    )


def read_gene_pathway_frame(file_gene_pathway: Path) -> pl.DataFrame:
    df = _read_two_col_tsv(
        file_gene_pathway,
        columns=["kegg_gene_id", "kegg_pathway_id"],
        schema={"kegg_gene_id": pl.String, "kegg_pathway_id": pl.String},
    )
    if df.height == 0:
        return df.with_columns(pl.lit(None, dtype=pl.String).alias("pathway_map_id"))
    return (
        df.with_columns(pl.col("kegg_pathway_id").str.strip_prefix("path:"))
        .with_columns(derive_pathway_map_id_expr().alias("pathway_map_id"))
        .unique()
        .sort("kegg_gene_id", "kegg_pathway_id")
    )


def read_gene_list_frame(file_gene_list: Path | None) -> pl.DataFrame:
    schema = {
        "kegg_gene_id": pl.String,
        "gene_symbol": pl.String,
        "gene_description": pl.String,
    }
    if file_gene_list is None:
        return pl.DataFrame(schema=schema)
    df = _read_two_col_tsv(
        file_gene_list,
        columns=["kegg_gene_id", "_gene_text"],
        schema={"kegg_gene_id": pl.String, "_gene_text": pl.String},
    )
    if df.height == 0:
        return pl.DataFrame(schema=schema)

    rows: list[dict[str, str | None]] = []
    for row in df.iter_rows(named=True):
        text = str(row["_gene_text"]).strip()
        if ";" in text:
            symbol, description = text.split(";", 1)
            symbol = symbol.strip() or None
            description = description.strip() or None
        else:
            symbol = text or None
            description = None
        rows.append(
            {
                "kegg_gene_id": row["kegg_gene_id"],
                "gene_symbol": symbol,
                "gene_description": description,
            }
        )
    return pl.DataFrame(rows, schema=schema).unique().sort("kegg_gene_id")


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
                df_conv_uniprot.select("kegg_gene_id"),
                df_conv_ncbi_geneid.select("kegg_gene_id"),
                df_gene_ko.select("kegg_gene_id"),
                df_gene_pathway.select("kegg_gene_id"),
                df_gene_list.select("kegg_gene_id"),
            ],
            how="vertical_relaxed",
        )
        .drop_nulls("kegg_gene_id")
        .unique()
    )
    if df_gene_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_MAPPING)

    validate_organism_code(df_gene_ids, organism_code=organism_code)
    return (
        df_gene_ids.join(df_conv_uniprot, on="kegg_gene_id", how="left")
        .join(df_conv_ncbi_geneid, on="kegg_gene_id", how="left")
        .join(df_gene_ko, on="kegg_gene_id", how="left")
        .join(df_gene_pathway, on="kegg_gene_id", how="left")
        .join(df_gene_list, on="kegg_gene_id", how="left")
        .with_columns(pl.lit(organism_code).alias("organism_code"))
        .select(COLS_MAPPING)
        .unique()
        .sort(COLS_MAPPING)
    )


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    namespace: KEGGNamespace,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    validate_namespace(namespace)
    col_join = column_by_namespace(namespace)
    cols_group = list(cols_group_id)
    cols_out = cols_group + ["input_id", "input_namespace"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="input_id",
            right_on=col_join,
            how="inner",
        )
        .with_columns(
            pl.col("input_id").alias(col_join),
            pl.lit(namespace).alias("input_namespace"),
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def extract_unmatched_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    cols_index = list(cols_group_id) + ["input_id"]
    df_mapped_input_ids = df_mapping.select(cols_index).unique().sort(cols_index)
    return (
        df_input_ids.join(df_mapped_input_ids, on=cols_index, how="anti")
        .select(cols_index)
        .sort(cols_index)
    )


def validate_namespace(namespace: str) -> None:
    if namespace not in NAMESPACE_VALUES:
        raise ValueError(
            "namespace must be one of: "
            f"{', '.join(NAMESPACE_VALUES)}; got {namespace!r}"
        )


def column_by_namespace(namespace: KEGGNamespace) -> str:
    match namespace:
        case "uniprot":
            return "uniprot_id"
        case "ncbi_gene":
            return "ncbi_gene_id"
        case "kegg_gene":
            return "kegg_gene_id"
        case _:
            validate_namespace(namespace)
            raise AssertionError("unreachable")


def validate_organism_code(df_gene_ids: pl.DataFrame, *, organism_code: str) -> None:
    prefix = f"{organism_code}:"
    df_invalid = df_gene_ids.filter(~pl.col("kegg_gene_id").str.starts_with(prefix))
    if df_invalid.height == 0:
        return
    examples = df_invalid.get_column("kegg_gene_id").head(5).to_list()
    raise ValueError(
        "KEGG gene IDs do not match organism_code: "
        f"organism_code={organism_code!r}, examples={examples!r}"
    )


def derive_pathway_map_id_expr() -> pl.Expr:
    return pl.lit("map") + pl.col("kegg_pathway_id").str.extract(
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
