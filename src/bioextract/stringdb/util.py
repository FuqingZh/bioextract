from collections.abc import Iterable, Mapping
from pathlib import Path

import polars as pl
import polars.selectors as pl_sel

from bioextract.stringdb.spec import GroupInputFrames, StringResourceLimits

from .constant import (
    SCHEMA_ALIASES,
    SCHEMA_EDGES,
    SCHEMA_GROUP_STRING_MAPPING,
    SCHEMA_GROUPS,
    SCHEMA_LINKS,
    SCHEMA_GROUP_EDGES,
    SCHEMA_GROUP_INPUT_IDS,
    RE_UNIPROT_PIPE,
    SCHEMA_PROTEIN_MAP,
    SCHEMA_UNMAPPED,
)


def _normalize_input_id(value: str) -> str:
    value = value.strip()
    if (match_pipe := RE_UNIPROT_PIPE.match(value)) is not None:
        return match_pipe.group(1).strip()
    return value


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


def validate_file_size(
    *,
    file_path: Path,
    size_max: int | None,
    label: str,
) -> None:
    if size_max is None:
        return
    file_size = file_path.stat().st_size
    if file_size > size_max:
        raise ValueError(
            f"{label} exceeds configured size limit: "
            f"path={file_path}, size_bytes={file_size}, limit_bytes={size_max}"
        )


def validate_input_id_count(
    *,
    num_input_ids: int,
    limits: StringResourceLimits,
) -> None:
    if limits.num_input_ids_max is None:
        return
    if num_input_ids > limits.num_input_ids_max:
        raise ValueError(
            "Normalized input ID count exceeds configured limit: "
            f"count={num_input_ids}, limit={limits.num_input_ids_max}"
        )


def validate_group_count(
    *,
    num_groups: int,
    limits: StringResourceLimits,
) -> None:
    if limits.num_groups_max is None:
        return
    if num_groups > limits.num_groups_max:
        raise ValueError(
            f"Group count exceeds configured limit: count={num_groups}, "
            f"limit={limits.num_groups_max}"
        )


def validate_group_ids(group_ids: list[str]) -> None:
    group_ids_seen: set[str] = set()
    for _id in group_ids:
        if not _id:
            raise ValueError("GroupId must be a non-empty string after normalization")
        if _id in group_ids_seen:
            raise ValueError(
                f"GroupId values must be unique after normalization: {_id!r}"
            )
        group_ids_seen.add(_id)


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


def validate_links_columns(cols_available: list[str], *, version: str) -> None:
    cols_required = set(SCHEMA_LINKS[version])
    if not cols_required.issubset(cols_available):
        raise ValueError(
            f"STRING links file for version {version} must contain columns "
            f"{sorted(cols_required)}; available={cols_available}"
        )


def create_input_id_frame(input_ids: Iterable[str]) -> pl.DataFrame:
    ids_normalized: list[str] = []
    for _id in input_ids:
        if input_id_normalized := _normalize_input_id(str(_id)):
            ids_normalized.append(input_id_normalized)

    if not ids_normalized:
        return pl.DataFrame(schema=SCHEMA_UNMAPPED)

    return (
        pl.DataFrame({"InputId": ids_normalized}, schema=SCHEMA_UNMAPPED)
        .unique(subset=["InputId"])
        .sort("InputId")
    )


def create_group_input_frames(
    group_to_ids: Mapping[str, Iterable[str]],
) -> GroupInputFrames:
    group_ids_normalized: list[str] = []
    group_ids_col: list[str] = []
    input_ids_col: list[str] = []

    for _grp_id, _ids in group_to_ids.items():
        group_id = str(_grp_id).strip()
        group_ids_normalized.append(group_id)
        for _id in _ids:
            if input_id_normalized := _normalize_input_id(str(_id)):
                group_ids_col.append(group_id)
                input_ids_col.append(input_id_normalized)

    validate_group_ids(group_ids_normalized)

    if group_ids_normalized:
        df_groups = pl.DataFrame(
            {"GroupId": group_ids_normalized}, schema=SCHEMA_GROUPS
        ).sort("GroupId")
    else:
        df_groups = pl.DataFrame(schema=SCHEMA_GROUPS)

    if not input_ids_col:
        return GroupInputFrames(
            df_groups=df_groups,
            df_input_ids=pl.DataFrame(schema=SCHEMA_GROUP_INPUT_IDS),
        )

    df_group_input_ids = (
        pl.DataFrame(
            {
                "GroupId": group_ids_col,
                "InputId": input_ids_col,
            },
            schema=SCHEMA_GROUP_INPUT_IDS,
        )
        .unique()
        .sort("GroupId", "InputId")
    )
    return GroupInputFrames(df_groups=df_groups, df_input_ids=df_group_input_ids)


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


def extract_group_string_mapping_frame(
    *,
    file_aliases: Path,
    df_input_ids: pl.DataFrame,
    source_rank_map: Mapping[str, int],
    col_string_id_aliases: str,
    has_source_aliases: bool,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_input_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_GROUP_STRING_MAPPING)

    lf_aliases = scan_aliases(file_aliases, version=version)
    lf_input_ids = df_input_ids.lazy()
    return (
        create_string_mapping_lazy_frame(
            lf_aliases=lf_aliases,
            lf_input_ids=lf_input_ids,
            source_rank_map=source_rank_map,
            col_string_id_aliases=col_string_id_aliases,
            has_source_aliases=has_source_aliases,
            cols_partition=["GroupId", "InputId", "StringId"],
            cols_sort_prefix=["GroupId", "InputId", "StringId"],
            cols_select_out=["GroupId", "InputId", "StringId", "MapSource"],
        )
        .sort(["GroupId", "InputId", "StringId", "MapSource"])
        .collect()
    )


def create_edges_lazy_frame(
    *,
    lf_links: pl.LazyFrame,
    lf_string_ids_a: pl.LazyFrame,
    lf_string_ids_b: pl.LazyFrame,
    thr_score_min: int,
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
            & pl.col("Score").ge(int(thr_score_min))
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


def extract_edges_frame(
    *,
    file_links: Path,
    df_string_ids: pl.DataFrame,
    thr_score_min: int,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_string_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_EDGES)

    lf_links = scan_links(file_links, version=version)
    cols_available = lf_links.collect_schema().names()
    validate_links_columns(cols_available, version=version)
    return (
        create_edges_lazy_frame(
            lf_links=lf_links,
            lf_string_ids_a=df_string_ids.lazy().rename({"StringId": "StringIdA"}),
            lf_string_ids_b=df_string_ids.lazy().rename({"StringId": "StringIdB"}),
            thr_score_min=thr_score_min,
            cols_join_left_a="StringIdA",
            cols_join_right_a="StringIdA",
            cols_join_left_b="StringIdB",
            cols_join_right_b="StringIdB",
            cols_partition=["_Lo", "_Hi"],
            cols_select_out=["StringIdA", "StringIdB", "Score"],
        )
        .sort(["StringIdA", "StringIdB"])
        .collect()
    )


def extract_group_edges_frame(
    *,
    file_links: Path,
    df_string_ids: pl.DataFrame,
    thr_score_min: int,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_string_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_GROUP_EDGES)

    lf_links = scan_links(file_links, version=version)
    cols_available = lf_links.collect_schema().names()
    validate_links_columns(cols_available, version=version)
    return (
        create_edges_lazy_frame(
            lf_links=lf_links,
            lf_string_ids_a=df_string_ids.lazy().rename({"StringId": "StringIdA"}),
            lf_string_ids_b=df_string_ids.lazy().rename({"StringId": "StringIdB"}),
            thr_score_min=thr_score_min,
            cols_join_left_a="StringIdA",
            cols_join_right_a="StringIdA",
            cols_join_left_b=["GroupId", "StringIdB"],
            cols_join_right_b=["GroupId", "StringIdB"],
            cols_partition=["GroupId", "_Lo", "_Hi"],
            cols_select_out=["GroupId", "StringIdA", "StringIdB", "Score"],
        )
        .sort(["GroupId", "StringIdA", "StringIdB"])
        .collect()
    )

