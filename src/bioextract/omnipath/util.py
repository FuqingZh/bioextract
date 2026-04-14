from pathlib import Path

import polars as pl
import polars.selectors as pl_sel

from bioextract._shared import validate_required_cols

from .constant import (
    COLS_RENAMED_ENZSUB,
    COLS_RENAMED_INTERACTIONS,
    SCHEMA_ENZSUB,
    SCHEMA_ENZSUB_RAW,
    SCHEMA_GROUP_ENZSUB,
    SCHEMA_GROUP_INTERACTIONS,
    SCHEMA_INTERACTIONS,
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
        lf_enzsub.select(
            pl_sel.by_name(COLS_RENAMED_ENZSUB)
            .cast(pl.String)
            .str.strip_chars()
            .replace("", None)
        )
        .with_columns(
            pl.col("residue_type").str.to_uppercase(),
            pl.col("modification").str.to_lowercase(),
        )
        .rename(COLS_RENAMED_ENZSUB)
        .drop_nulls(["SourceId", "TargetId", "_ResidueType", "_ResidueOffset"])
        .with_columns(
            pl.concat_str(["_ResidueType", "_ResidueOffset"]).alias("TargetSite")
        )
        .drop(["_ResidueType", "_ResidueOffset"])
    )


def _create_interactions_base_lazy_frame(lf_interactions: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf_interactions.select(
            [
                pl_sel.by_name("source", "target")
                .cast(pl.String)
                .str.strip_chars()
                .replace("", None)
                .name.keep(),
                _normalize_bool_expr(
                    "is_directed",
                    "is_stimulation",
                    "is_inhibition",
                ),
            ]
        )
        .rename(COLS_RENAMED_INTERACTIONS)
        .drop_nulls(["SourceId", "TargetId"])
    )


def _create_validated_enzsub_base_lazy_frame(file_enzsub: Path) -> pl.LazyFrame:
    lf_enzsub = scan_enzsub(file_enzsub)
    validate_required_cols(
        cols_available=lf_enzsub.collect_schema().names(),
        cols_required=SCHEMA_ENZSUB_RAW.keys(),
        context="OmniPath enzsub file",
    )
    return _create_enzsub_base_lazy_frame(lf_enzsub)


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


def _has_any_rows(lf: pl.LazyFrame) -> bool:
    return lf.limit(1).collect().height > 0


def has_any_enzsub_relation(*, file_enzsub: Path) -> bool:
    return _has_any_rows(_create_validated_enzsub_base_lazy_frame(file_enzsub))


def has_any_interaction_relation(*, file_interactions: Path) -> bool:
    return _has_any_rows(
        _create_validated_interactions_base_lazy_frame(file_interactions)
    )


def has_any_enzsub_modification(*, file_enzsub: Path, modification: str) -> bool:
    modification_normalized = str(modification).strip().lower()
    if not modification_normalized:
        raise ValueError("Modification must be a non-empty string after normalization")

    return _has_any_rows(
        _create_validated_enzsub_base_lazy_frame(file_enzsub).filter(
            pl.col("Modification").eq(modification_normalized)
        )
    )


def _concat_matched_frames(
    *,
    lf_left: pl.LazyFrame,
    lf_right: pl.LazyFrame,
    cols_unique: list[str],
    cols_sort: list[str],
    cols_select: list[str],
) -> pl.DataFrame:
    return (
        pl.concat([lf_left, lf_right], how="vertical_relaxed")
        .unique(subset=cols_unique)
        .select(cols_select)
        .sort(cols_sort)
        .collect()
    )


def extract_enzsub_frame(
    *,
    file_enzsub: Path,
    df_input_ids: pl.DataFrame,
    cols_group_id: tuple[str, ...] = (),
) -> pl.DataFrame:
    if df_input_ids.height == 0:
        return pl.DataFrame(
            schema=SCHEMA_GROUP_ENZSUB if cols_group_id else SCHEMA_ENZSUB
        )

    lf_enzsub = scan_enzsub(file_enzsub)
    validate_required_cols(
        cols_available=lf_enzsub.collect_schema().names(),
        cols_required=SCHEMA_ENZSUB_RAW.keys(),
        context="OmniPath enzsub file",
    )
    lf_base = _create_enzsub_base_lazy_frame(lf_enzsub)

    lf_source = lf_base.join(
        df_input_ids.lazy().rename({"InputId": "SourceId"}),
        on="SourceId",
        how="inner",
    )
    lf_target = lf_base.join(
        df_input_ids.lazy().rename({"InputId": "TargetId"}),
        on="TargetId",
        how="inner",
    )
    cols_group = list(cols_group_id)
    return _concat_matched_frames(
        lf_left=lf_source,
        lf_right=lf_target,
        cols_unique=cols_group + ["SourceId", "TargetId", "TargetSite", "Modification"],
        cols_sort=cols_group + ["SourceId", "TargetId", "TargetSite", "Modification"],
        cols_select=cols_group + ["SourceId", "TargetId", "TargetSite", "Modification"],
    )


def extract_interactions_frame(
    *,
    file_interactions: Path,
    df_input_ids: pl.DataFrame,
    cols_group_id: tuple[str, ...] = (),
) -> pl.DataFrame:
    if df_input_ids.height == 0:
        return pl.DataFrame(
            schema=SCHEMA_GROUP_INTERACTIONS if cols_group_id else SCHEMA_INTERACTIONS
        )

    lf_interactions = scan_interactions(file_interactions)
    validate_required_cols(
        cols_available=lf_interactions.collect_schema().names(),
        cols_required=SCHEMA_INTERACTIONS_RAW.keys(),
        context="OmniPath interactions file",
    )
    lf_base = _create_interactions_base_lazy_frame(lf_interactions)

    lf_source = lf_base.join(
        df_input_ids.lazy().rename({"InputId": "SourceId"}),
        on="SourceId",
        how="inner",
    )
    lf_target = lf_base.join(
        df_input_ids.lazy().rename({"InputId": "TargetId"}),
        on="TargetId",
        how="inner",
    )
    cols_group = list(cols_group_id)
    return _concat_matched_frames(
        lf_left=lf_source,
        lf_right=lf_target,
        cols_unique=cols_group + [
            "SourceId",
            "TargetId",
            "IsDirected",
            "IsStimulation",
            "IsInhibition",
        ],
        cols_sort=cols_group + [
            "SourceId",
            "TargetId",
            "IsDirected",
            "IsStimulation",
            "IsInhibition",
        ],
        cols_select=cols_group + [
            "SourceId",
            "TargetId",
            "IsDirected",
            "IsStimulation",
            "IsInhibition",
        ],
    )
