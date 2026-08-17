from pathlib import Path

import polars as pl
import polars.selectors as pl_sel

from bioextract._shared import validate_required_cols

from .constant import (
    SCHEMA_ENZSUB_RAW,
    SCHEMA_INTERACTIONS_RAW,
)


def _normalize_bool_expr(*col_names: str) -> pl.Expr:
    expr_raw = (
        pl_sel.by_name(*col_names).cast(pl.String).str.strip_chars().str.to_lowercase()
    )
    return (
        pl.when(expr_raw.is_in(["true", "1"]))
        .then(pl.lit(True))
        .when(expr_raw.is_in(["false", "0"]))
        .then(pl.lit(False))
        .otherwise(None)
        .cast(pl.Boolean)
        .name.keep()
    )


def scan_enzsub(file_enzsub: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_enzsub,
        separator="\t",
        schema_overrides=SCHEMA_ENZSUB_RAW,
        null_values=["", "-"],
    )


def scan_interactions(file_interactions: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_interactions,
        separator="\t",
        schema_overrides=SCHEMA_INTERACTIONS_RAW,
        null_values=["", "-"],
    )


def _create_enzsub_base_lazy_frame(lf_enzsub: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf_enzsub.select(pl_sel.by_name(SCHEMA_ENZSUB_RAW).cast(pl.String).name.keep())
        .with_columns(
            pl.col("enzyme").str.strip_chars().replace("", None).alias("_enzyme"),
            pl.col("substrate").str.strip_chars().replace("", None).alias("_substrate"),
            pl.concat_str(
                [
                    pl.col("residue_type").str.strip_chars(),
                    pl.col("residue_offset").str.strip_chars(),
                ],
                separator="",
                ignore_nulls=True,
            )
            .replace("", None)
            .alias("target_site"),
        )
        .drop_nulls(["_enzyme", "_substrate"])
    )


def _create_interactions_base_lazy_frame(lf_interactions: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf_interactions.select(
            [
                pl_sel.by_name("source", "target").cast(pl.String).name.keep(),
                _normalize_bool_expr(
                    "is_directed",
                    "is_stimulation",
                    "is_inhibition",
                ),
            ]
        )
        .with_columns(
            pl.col("source").str.strip_chars().replace("", None).alias("_source"),
            pl.col("target").str.strip_chars().replace("", None).alias("_target"),
        )
        .drop_nulls(["_source", "_target"])
    )


def _create_validated_enzsub_base_lazy_frame(file_enzsub: Path) -> pl.LazyFrame:
    lf_enzsub = scan_enzsub(file_enzsub)
    validate_required_cols(
        cols_available=lf_enzsub.collect_schema().names(),
        cols_required=SCHEMA_ENZSUB_RAW.keys(),
        context="OmniPath enzsub file",
    )
    return _create_enzsub_base_lazy_frame(lf_enzsub)


def create_enzsub_lazy_frame(file_enzsub: Path) -> pl.LazyFrame:
    """Return the validated, normalized enzyme-substrate scan lazily."""
    return _create_validated_enzsub_base_lazy_frame(file_enzsub)


def _create_validated_interactions_base_lazy_frame(
    file_interactions: Path,
) -> pl.LazyFrame:
    lf_interactions = scan_interactions(file_interactions)
    validate_required_cols(
        cols_available=lf_interactions.collect_schema().names(),
        cols_required=SCHEMA_INTERACTIONS_RAW.keys(),
        context="OmniPath interactions file",
    )
    return _create_interactions_base_lazy_frame(lf_interactions)


def create_interactions_lazy_frame(file_interactions: Path) -> pl.LazyFrame:
    """Return the validated, normalized interaction scan lazily."""
    return _create_validated_interactions_base_lazy_frame(file_interactions)
