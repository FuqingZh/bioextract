from collections.abc import Mapping
from pathlib import Path

import polars as pl
import polars.selectors as pl_sel

from bioextract._shared import (
    RE_UNIPROT_PIPE,
)

from .constant import (
    SCHEMA_ALIASES,
    SCHEMA_LINKS,
    SCHEMA_PROTEIN_MAP,
)


def _normalize_input_id_expr(col_name: str) -> pl.Expr:
    expr_raw = pl.col(col_name).cast(pl.String).str.strip_chars()
    return (
        pl.coalesce(
            expr_raw.str.extract(RE_UNIPROT_PIPE.pattern, 1),
            expr_raw,
        )
        .replace("", None)
        .alias("InputId")
    )


def scan_aliases(file_aliases: Path, version: str) -> pl.LazyFrame:
    return pl.scan_csv(
        file_aliases,
        separator="\t",
        schema_overrides=SCHEMA_ALIASES[version],
        null_values=["", "-"],
    )


def scan_links(file_links: Path, version: str) -> pl.LazyFrame:
    return pl.scan_csv(
        file_links,
        separator=" ",
        schema_overrides=SCHEMA_LINKS[version],
        null_values=["", "-"],
    )


def validate_alias_required_cols(cols_available: list[str], version: str) -> None:
    cols_expected = set(SCHEMA_ALIASES[version])
    if "alias" not in cols_expected or "alias" not in cols_available:
        raise ValueError(
            f"STRING aliases file for version {version} must contain 'alias'; "
            f"available={cols_available}"
        )


def infer_alias_id_col(cols_available: list[str], *, version: str) -> str:
    cols_expected = set(SCHEMA_ALIASES[version])

    if "#string_protein_id" in cols_expected and "#string_protein_id" in cols_available:
        col_string_id = "#string_protein_id"
    elif "string_protein_id" in cols_expected and "string_protein_id" in cols_available:
        col_string_id = "string_protein_id"
    else:
        raise ValueError(
            f"STRING aliases file for version {version} must contain "
            f"'#string_protein_id' or 'string_protein_id'; available={cols_available}"
        )

    return col_string_id


def create_unknown_source_rank(base: int, add: int = 100) -> int:
    return base + add


def create_string_mapping_lazy_frame(
    *,
    lf_aliases: pl.LazyFrame,
    lf_input_ids: pl.LazyFrame,
    source_rank_map: Mapping[str, int],
    col_string_id_aliases: str,
    has_source_aliases: bool,
    cols_partition: list[str],
    cols_sort_prefix: list[str],
    cols_select_out: list[str],
) -> pl.LazyFrame:
    lf_string_mapping = (
        lf_aliases.select(
            [
                pl.col(col_string_id_aliases)
                .cast(pl.String)
                .str.strip_chars()
                .replace("", None)
                .alias("StringId"),
                _normalize_input_id_expr("alias"),
                (
                    pl.col("source").cast(pl.String).str.strip_chars()
                    if has_source_aliases
                    else pl.lit("unknown")
                ).alias("MapSource"),
            ]
        )
        .drop_nulls(["InputId", "StringId"])
        .join(lf_input_ids, on="InputId", how="inner")
    )

    rank_default = create_unknown_source_rank(len(source_rank_map))
    if source_rank_map:
        lf_rank = pl.LazyFrame(
            {
                "MapSource": list(source_rank_map.keys()),
                "_SourceRank": list(source_rank_map.values()),
            },
            schema={"MapSource": pl.String, "_SourceRank": pl.Int64},
        )
        lf_string_mapping = lf_string_mapping.join(lf_rank, on="MapSource", how="left")
    else:
        lf_string_mapping = lf_string_mapping.with_columns(
            pl.lit(None).cast(pl.Int64).alias("_SourceRank")
        )

    return (
        lf_string_mapping.with_columns(
            pl.col("_SourceRank").fill_null(rank_default).cast(pl.Int64)
        )
        .sort([*cols_sort_prefix, "_SourceRank", "MapSource"])
        .unique(subset=cols_partition, keep="first")
        .drop("_SourceRank")
        .select(cols_select_out)
    )


def extract_string_mapping_frame(
    *,
    file_aliases: Path,
    df_input_ids: pl.DataFrame,
    source_rank_map: Mapping[str, int],
    col_string_id_aliases: str,
    has_source_aliases: bool,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_input_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_PROTEIN_MAP)

    lf_aliases = scan_aliases(file_aliases, version=version)
    lf_input_ids = df_input_ids.lazy()
    return (
        create_string_mapping_lazy_frame(
            lf_aliases=lf_aliases,
            lf_input_ids=lf_input_ids,
            source_rank_map=source_rank_map,
            col_string_id_aliases=col_string_id_aliases,
            has_source_aliases=has_source_aliases,
            cols_partition=["InputId", "StringId"],
            cols_sort_prefix=["InputId", "StringId"],
            cols_select_out=["InputId", "StringId", "MapSource"],
        )
        .sort(["InputId", "StringId", "MapSource"])
        .collect()
    )


def create_edges_lazy_frame(
    *,
    lf_links: pl.LazyFrame,
    lf_string_ids_a: pl.LazyFrame,
    lf_string_ids_b: pl.LazyFrame,
    min_combined_score: int,
    cols_join_left_a: str | list[str],
    cols_join_right_a: str | list[str],
    cols_join_left_b: str | list[str],
    cols_join_right_b: str | list[str],
    cols_partition: list[str],
    cols_select_out: list[str],
) -> pl.LazyFrame:
    return (
        lf_links.select(
            pl.col("protein1").cast(pl.String).str.strip_chars().alias("StringIdA"),
            pl.col("protein2").cast(pl.String).str.strip_chars().alias("StringIdB"),
            pl.col("combined_score").cast(pl.Int64).alias("Score"),
        )
        .with_columns(pl_sel.by_name("StringIdA", "StringIdB").replace("", None))
        .drop_nulls(["StringIdA", "StringIdB", "Score"])
        .filter(
            pl.col("StringIdA").ne(pl.col("StringIdB"))
            & pl.col("Score").ge(int(min_combined_score))
        )
        .join(
            lf_string_ids_a,
            left_on=cols_join_left_a,
            right_on=cols_join_right_a,
            how="inner",
        )
        .join(
            lf_string_ids_b,
            left_on=cols_join_left_b,
            right_on=cols_join_right_b,
            how="inner",
        )
        .with_columns(
            pl.when(pl.col("StringIdA").le(pl.col("StringIdB")))
            .then(pl.col("StringIdA"))
            .otherwise(pl.col("StringIdB"))
            .alias("_Lo"),
            pl.when(pl.col("StringIdA").le(pl.col("StringIdB")))
            .then(pl.col("StringIdB"))
            .otherwise(pl.col("StringIdA"))
            .alias("_Hi"),
        )
        .group_by(cols_partition)
        .agg(pl.col("Score").max().cast(pl.Int64).alias("Score"))
        .rename({"_Lo": "StringIdA", "_Hi": "StringIdB"})
        .select(cols_select_out)
    )
