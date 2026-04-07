from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import polars as pl
import polars.selectors as pl_sel

__all__ = [
    "StringDb",
    "StringResourceLimits",
]

_DEFAULT_SOURCE_RANK_MAP: dict[str, int] = {
    "UniProt_AC": 1,
    "UniProt_ID": 2,
    "UniProt_GN_Name": 3,
    "UniProt_GN_Synonyms": 4,
    "UniProt_DR_GeneID": 5,
    "KEGG_GENEID": 6,
}

SCHEMA_ALIASES = {
    "v12.0": {
        "#string_protein_id": pl.String,
        "string_protein_id": pl.String,
        "alias": pl.String,
        "source": pl.String,
    },
}
SCHEMA_LINKS = {
    "v12.0": {
        "protein1": pl.String,
        "protein2": pl.String,
        "combined_score": pl.Int64,
    },
}
SCHEMA_PROTEIN_MAP = {
    "InputId": pl.String,
    "StringId": pl.String,
    "MapSource": pl.String,
}
SCHEMA_EDGES = {
    "StringIdA": pl.String,
    "StringIdB": pl.String,
    "Score": pl.Int64,
}
SCHEMA_UNMAPPED = {
    "InputId": pl.String,
}
SCHEMA_GROUPS = {
    "GroupId": pl.String,
}
SCHEMA_GROUP_INPUT_IDS = {
    "GroupId": pl.String,
    "InputId": pl.String,
}
SCHEMA_GROUP_STRING_MAPPING = {
    "GroupId": pl.String,
    "InputId": pl.String,
    "StringId": pl.String,
    "MapSource": pl.String,
}
SCHEMA_GROUP_UNMAPPED = {
    "GroupId": pl.String,
    "InputId": pl.String,
}
SCHEMA_GROUP_STRING_IDS = {
    "GroupId": pl.String,
    "StringId": pl.String,
}
SCHEMA_GROUP_EDGES = {
    "GroupId": pl.String,
    "StringIdA": pl.String,
    "StringIdB": pl.String,
    "Score": pl.Int64,
}
SCHEMA_GROUP_METRICS = {
    "GroupId": pl.String,
    "NumInputIds": pl.Int64,
    "NumMappedIds": pl.Int64,
    "NumUnmappedIds": pl.Int64,
    "NumStringIds": pl.Int64,
    "NumEdges": pl.Int64,
}

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")


@dataclass(frozen=True, slots=True)
class StringResourceLimits:
    file_aliases_bytes_max: int | None = 512 * 1024 * 1024
    file_links_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000


@dataclass(frozen=True, slots=True)
class StringExtractMetrics:
    num_input_ids: int
    num_mapped_ids: int
    num_unmapped_ids: int
    num_string_ids: int
    num_edges: int


@dataclass(frozen=True, slots=True)
class StringExtractResult:
    df_edges: pl.DataFrame
    df_protein_map: pl.DataFrame
    df_unmapped: pl.DataFrame
    metrics: StringExtractMetrics


@dataclass(frozen=True, slots=True)
class _StringAliasSchema:
    col_string_id: str
    has_source: bool


def _normalize_input_id(value: str) -> str:
    value = value.strip()
    match_pipe = RE_UNIPROT_PIPE.match(value)
    if match_pipe is not None:
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


def _scan_aliases(file_aliases: Path, version: str) -> pl.LazyFrame:
    return pl.scan_csv(
        file_aliases,
        separator="\t",
        schema_overrides=SCHEMA_ALIASES[version],
        null_values=["", "-"],
    )


def _scan_links(file_links: Path, version: str) -> pl.LazyFrame:
    return pl.scan_csv(
        file_links,
        separator=" ",
        schema_overrides=SCHEMA_LINKS[version],
        null_values=["", "-"],
    )


def _validate_file_size(
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


def _validate_input_id_count(
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


def _validate_group_count(
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


def _validate_group_ids(group_ids: list[str]) -> None:
    group_ids_seen: set[str] = set()
    for group_id in group_ids:
        if not group_id:
            raise ValueError("GroupId must be a non-empty string after normalization")
        if group_id in group_ids_seen:
            raise ValueError(
                f"GroupId values must be unique after normalization: {group_id!r}"
            )
        group_ids_seen.add(group_id)


def _validate_alias_columns(
    cols_available: list[str], *, version: str
) -> _StringAliasSchema:
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

    if "alias" not in cols_expected or "alias" not in cols_available:
        raise ValueError(
            f"STRING aliases file for version {version} must contain 'alias'; "
            f"available={cols_available}"
        )
    return _StringAliasSchema(
        col_string_id=col_string_id,
        has_source="source" in cols_available,
    )


def _validate_links_columns(cols_available: list[str], *, version: str) -> None:
    cols_required = set(SCHEMA_LINKS[version])
    if not cols_required.issubset(cols_available):
        raise ValueError(
            f"STRING links file for version {version} must contain columns "
            f"{sorted(cols_required)}; available={cols_available}"
        )


def _create_input_id_frame(input_ids: Iterable[str]) -> pl.DataFrame:
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


def _create_group_input_frames(
    group_to_ids: Mapping[str, Iterable[str]],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    group_ids_normalized: list[str] = []
    rows_group_input_ids: list[dict[str, str]] = []

    for group_id_raw, input_ids in group_to_ids.items():
        group_id = str(group_id_raw).strip()
        group_ids_normalized.append(group_id)
        for _id in input_ids:
            if input_id_normalized := _normalize_input_id(str(_id)):
                rows_group_input_ids.append(
                    {"GroupId": group_id, "InputId": input_id_normalized}
                )

    _validate_group_ids(group_ids_normalized)

    if group_ids_normalized:
        df_groups = (
            pl.DataFrame({"GroupId": group_ids_normalized}, schema=SCHEMA_GROUPS)
            .unique(subset=["GroupId"])
            .sort("GroupId")
        )
    else:
        df_groups = pl.DataFrame(schema=SCHEMA_GROUPS)

    if not rows_group_input_ids:
        return df_groups, pl.DataFrame(schema=SCHEMA_GROUP_INPUT_IDS)

    df_group_input_ids = (
        pl.DataFrame(rows_group_input_ids, schema=SCHEMA_GROUP_INPUT_IDS)
        .unique(subset=["GroupId", "InputId"])
        .sort(["GroupId", "InputId"])
    )
    return df_groups, df_group_input_ids


def _create_unknown_source_rank(base: int, add: int = 100) -> int:
    return base + add


def _create_string_mapping_lazy_frame(
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

    rank_default = _create_unknown_source_rank(len(source_rank_map))
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


def _extract_string_mapping_frame(
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

    lf_aliases = _scan_aliases(file_aliases, version=version)
    lf_input_ids = df_input_ids.lazy()
    return (
        _create_string_mapping_lazy_frame(
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


def _extract_group_string_mapping_frame(
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

    lf_aliases = _scan_aliases(file_aliases, version=version)
    lf_input_ids = df_input_ids.lazy()
    return (
        _create_string_mapping_lazy_frame(
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


def _create_edges_lazy_frame(
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


def _extract_edges_frame(
    *,
    file_links: Path,
    df_string_ids: pl.DataFrame,
    thr_score_min: int,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_string_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_EDGES)

    lf_links = _scan_links(file_links, version=version)
    cols_available = lf_links.collect_schema().names()
    _validate_links_columns(cols_available, version=version)
    return (
        _create_edges_lazy_frame(
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


def _extract_group_edges_frame(
    *,
    file_links: Path,
    df_string_ids: pl.DataFrame,
    thr_score_min: int,
    version: str = "v12.0",
) -> pl.DataFrame:
    if df_string_ids.height == 0:
        return pl.DataFrame(schema=SCHEMA_GROUP_EDGES)

    lf_links = _scan_links(file_links, version=version)
    cols_available = lf_links.collect_schema().names()
    _validate_links_columns(cols_available, version=version)
    return (
        _create_edges_lazy_frame(
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


def _extract_group_metrics_frame(
    *,
    df_groups: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    df_string_mapping: pl.DataFrame,
    df_unmapped_input_ids: pl.DataFrame,
    df_string_ids: pl.DataFrame,
    df_edges: pl.DataFrame,
) -> pl.DataFrame:
    if df_groups.height == 0:
        return pl.DataFrame(schema=SCHEMA_GROUP_METRICS)

    df_num_input_ids = df_input_ids.group_by("GroupId").agg(
        pl.len().cast(pl.Int64).alias("NumInputIds")
    )
    df_num_mapped_ids = df_string_mapping.group_by("GroupId").agg(
        pl.col("InputId").n_unique().cast(pl.Int64).alias("NumMappedIds")
    )
    df_num_unmapped_ids = df_unmapped_input_ids.group_by("GroupId").agg(
        pl.len().cast(pl.Int64).alias("NumUnmappedIds")
    )
    df_num_string_ids = df_string_ids.group_by("GroupId").agg(
        pl.len().cast(pl.Int64).alias("NumStringIds")
    )
    df_num_edges = df_edges.group_by("GroupId").agg(
        pl.len().cast(pl.Int64).alias("NumEdges")
    )

    return (
        df_groups.join(df_num_input_ids, on="GroupId", how="left")
        .join(df_num_mapped_ids, on="GroupId", how="left")
        .join(df_num_unmapped_ids, on="GroupId", how="left")
        .join(df_num_string_ids, on="GroupId", how="left")
        .join(df_num_edges, on="GroupId", how="left")
        .with_columns(
            pl_sel.numeric().fill_null(0).cast(pl.Int64),
        )
        .select(
            [
                "GroupId",
                "NumInputIds",
                "NumMappedIds",
                "NumUnmappedIds",
                "NumStringIds",
                "NumEdges",
            ]
        )
        .sort("GroupId")
    )


@dataclass(slots=True)
class StringDb:
    """Path-first access to a local STRING snapshot.

    `StringDb` is the public entrypoint for extracting STRING mappings and
    interaction edges from local `protein.aliases` and `protein.links` files.
    It keeps dataset-level configuration such as source-rank policy, supported
    STRING version, and fail-fast resource limits.

    Construct instances with :meth:`from_files` and then choose one of two
    query styles:

    - single selection via :meth:`select_ids`
    - grouped selection via :meth:`select_groups`

    Both query styles preserve lazy scan behavior for the large STRING inputs
    until the final materialized result is requested.

    Examples:
        Single query:

            db = StringDb.from_files(
                file_aliases="9606.protein.aliases.v12.0.txt.gz",
                file_links="9606.protein.links.v12.0.txt.gz",
            )
            result = (
                db.select_ids(["TP53", "EGFR"])
                .with_score_min(400)
                .extract_result()
            )

        Grouped query:

            db = StringDb.from_files(
                file_aliases="9606.protein.aliases.v12.0.txt.gz",
                file_links="9606.protein.links.v12.0.txt.gz",
            )
            df_edges = (
                db.select_groups({"TumorA": ["TP53", "EGFR"], "TumorB": ["CDK2"]})
                .with_score_min(400)
                .extract_edges()
            )
    """

    file_aliases: Path
    file_links: Path
    source_rank_map: Mapping[str, int]
    limits: StringResourceLimits = field(default_factory=StringResourceLimits)
    version: Literal["v12.0"] = "v12.0"
    _alias_schema_cached: _StringAliasSchema | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_SOURCE_RANK_MAP = _DEFAULT_SOURCE_RANK_MAP
    DEFAULT_RESOURCE_LIMITS = StringResourceLimits()

    @classmethod
    def from_files(
        cls,
        file_aliases: str | Path,
        file_links: str | Path,
        *,
        source_rank_map: Mapping[str, int] = _DEFAULT_SOURCE_RANK_MAP,
        limits: StringResourceLimits | None = None,
        version: Literal["v12.0"] = "v12.0",
    ) -> StringDb:
        """Create a dataset handle from local STRING aliases and links files.

        Args:
            file_aliases: Path to a local STRING `protein.aliases` text or gzip
                file.
            file_links: Path to a local STRING `protein.links` text or gzip
                file.
            source_rank_map: Source-priority mapping used to break ties when the
                same input ID maps to the same STRING ID through multiple alias
                sources.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used.
            version: Declared STRING schema version used for CSV schema
                overrides and column validation.

        Returns:
            A dataset handle that can produce single or grouped selections.

        Raises:
            FileNotFoundError: If either input file does not exist.
            ValueError: If the requested version is unsupported or a configured
                file-size limit is exceeded.
        """
        file_aliases_path = Path(file_aliases)
        file_links_path = Path(file_links)
        if not file_aliases_path.exists():
            raise FileNotFoundError(
                f"STRING aliases file not found: {file_aliases_path}"
            )
        if not file_links_path.exists():
            raise FileNotFoundError(f"STRING links file not found: {file_links_path}")
        if version not in SCHEMA_ALIASES or version not in SCHEMA_LINKS:
            raise ValueError(f"Unsupported STRING version: {version}")

        limits_resolved = StringResourceLimits() if limits is None else limits
        _validate_file_size(
            file_path=file_aliases_path,
            size_max=limits_resolved.file_aliases_bytes_max,
            label="STRING aliases file",
        )
        _validate_file_size(
            file_path=file_links_path,
            size_max=limits_resolved.file_links_bytes_max,
            label="STRING links file",
        )

        return cls(
            file_aliases=file_aliases_path,
            file_links=file_links_path,
            source_rank_map=dict(source_rank_map),
            limits=limits_resolved,
            version=version,
        )

    def select_ids(self, input_ids: Iterable[str]) -> StringSelection:
        """Create a single-query selection from input IDs.

        Input IDs are normalized before selection:

        - surrounding whitespace is stripped
        - UniProt pipe-style values such as `sp|P04637|P53_HUMAN` are reduced
          to the middle accession token
        - empty normalized IDs are dropped
        - duplicates are removed

        Args:
            input_ids: Input protein, gene, or alias identifiers to resolve
                against the STRING aliases table.

        Returns:
            A cached selection object that can extract mapping, unmapped IDs,
            edges, or a bundled result.

        Raises:
            ValueError: If the normalized input-ID count exceeds the configured
                dataset limits.
        """
        df_input_ids = _create_input_id_frame(input_ids)
        _validate_input_id_count(
            num_input_ids=df_input_ids.height,
            limits=self.limits,
        )
        return StringSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            thr_score_min=0,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
    ) -> StringGroupSelection:
        """Create a grouped selection from multiple input-ID sets.

        Each group key is normalized with `str(...).strip()`. Input IDs within
        each group follow the same normalization rules as :meth:`select_ids`.
        Grouped extraction keeps groups isolated in the returned flat tables by
        carrying a `GroupId` column through mapping, unmapped, edge, and metric
        outputs.

        Args:
            group_to_ids: Mapping from group label to iterable of input
                identifiers.

        Returns:
            A grouped selection object that can extract grouped flat tables
            with `GroupId` as the leading column.

        Raises:
            ValueError: If any normalized `GroupId` is empty, if normalized
                `GroupId` values are duplicated, if the number of groups
                exceeds the configured limit, or if the total normalized
                input-ID count exceeds the configured limit.
        """
        df_groups, df_input_ids = _create_group_input_frames(group_to_ids)
        _validate_group_count(
            num_groups=df_groups.height,
            limits=self.limits,
        )
        _validate_input_id_count(
            num_input_ids=df_input_ids.height,
            limits=self.limits,
        )
        return StringGroupSelection(
            dataset=self,
            _df_groups=df_groups,
            _df_input_ids=df_input_ids,
            thr_score_min=0,
        )

    @property
    def alias_schema(self) -> _StringAliasSchema:
        if self._alias_schema_cached is not None:
            return self._alias_schema_cached
        lf_aliases = _scan_aliases(self.file_aliases, version=self.version)
        self._alias_schema_cached = _validate_alias_columns(
            lf_aliases.collect_schema().names(),
            version=self.version,
        )
        return self._alias_schema_cached


@dataclass(slots=True)
class StringSelection:
    dataset: StringDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    thr_score_min: int = 0
    _df_protein_map: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)
    _df_string_ids: pl.DataFrame | None = field(default=None, repr=False)
    _df_edges: pl.DataFrame | None = field(default=None, repr=False)
    _result: StringExtractResult | None = field(default=None, repr=False)

    def with_score_min(self, thr_score_min: int) -> StringSelection:
        """Create a new selection with a different minimum STRING score.

        Cached mapping-related frames are reused. Edge-related caches are not
        reused because the score threshold changes the edge result.

        Args:
            thr_score_min: Minimum `combined_score` required for retained
                STRING edges.

        Returns:
            A new selection sharing cached mapping state with the current
            selection.
        """
        return StringSelection(
            dataset=self.dataset,
            _df_input_ids=self._df_input_ids,
            thr_score_min=int(thr_score_min),
            _df_protein_map=self._df_protein_map,
            _df_unmapped=self._df_unmapped,
            _df_string_ids=self._df_string_ids,
        )

    def extract_string_mapping(self) -> pl.DataFrame:
        """Extract the input-to-STRING mapping table for this selection.

        Returns:
            A materialized table with columns `InputId`, `StringId`, and
            `MapSource`.

        Raises:
            ValueError: If the aliases file is missing required columns for the
                configured STRING version.
        """
        if self._df_protein_map is None:
            alias_schema = self.dataset.alias_schema
            self._df_protein_map = _extract_string_mapping_frame(
                file_aliases=self.dataset.file_aliases,
                df_input_ids=self._df_input_ids,
                source_rank_map=self.dataset.source_rank_map,
                col_string_id_aliases=alias_schema.col_string_id,
                has_source_aliases=alias_schema.has_source,
                version=self.dataset.version,
            )
        return self._df_protein_map

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that were not mapped to STRING IDs.

        Returns:
            A materialized table with one `InputId` column containing
            normalized inputs that did not resolve through the aliases table.
        """
        if self._df_unmapped is None:
            df_mapped_input_ids = self.extract_string_mapping().select("InputId")
            self._df_unmapped = (
                self._df_input_ids.join(
                    df_mapped_input_ids.unique(),
                    on="InputId",
                    how="anti",
                )
                .select(["InputId"])
                .sort("InputId")
            )
        return self._df_unmapped

    def extract_edges(self) -> pl.DataFrame:
        """Extract the STRING subnetwork induced by the mapped STRING IDs.

        Returns:
            A materialized edge table with columns `StringIdA`, `StringIdB`,
            and `Score`.

        Raises:
            ValueError: If the links file is missing required columns for the
                configured STRING version.
        """
        if self._df_edges is None:
            self._df_edges = _extract_edges_frame(
                file_links=self.dataset.file_links,
                df_string_ids=self._extract_string_ids(),
                thr_score_min=self.thr_score_min,
                version=self.dataset.version,
            )
        return self._df_edges

    def extract_result(self) -> StringExtractResult:
        """Extract the full single-query result bundle.

        The bundle contains the edge table, the input-to-STRING mapping table,
        the unmapped input table, and summary metrics derived from those
        materialized outputs.

        Returns:
            A bundled single-query extraction result with cached DataFrame
            members and summary counts.
        """
        if self._result is None:
            df_protein_map = self.extract_string_mapping()
            df_unmapped = self.extract_unmapped_input_ids()
            df_edges = self.extract_edges()
            df_string_ids = self._extract_string_ids()
            self._result = StringExtractResult(
                df_edges=df_edges,
                df_protein_map=df_protein_map,
                df_unmapped=df_unmapped,
                metrics=StringExtractMetrics(
                    num_input_ids=self._df_input_ids.height,
                    num_mapped_ids=df_protein_map.select("InputId").n_unique(),
                    num_unmapped_ids=df_unmapped.height,
                    num_string_ids=df_string_ids.height,
                    num_edges=df_edges.height,
                ),
            )
        return self._result

    def _extract_string_ids(self) -> pl.DataFrame:
        if self._df_string_ids is None:
            self._df_string_ids = (
                self.extract_string_mapping()
                .select("StringId")
                .unique()
                .sort("StringId")
            )
        return self._df_string_ids


@dataclass(slots=True)
class StringGroupSelection:
    dataset: StringDb
    _df_groups: pl.DataFrame = field(repr=False)
    _df_input_ids: pl.DataFrame = field(repr=False)
    thr_score_min: int = 0
    _df_string_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped_input_ids: pl.DataFrame | None = field(default=None, repr=False)
    _df_string_ids: pl.DataFrame | None = field(default=None, repr=False)
    _df_edges: pl.DataFrame | None = field(default=None, repr=False)
    _df_metrics: pl.DataFrame | None = field(default=None, repr=False)

    def with_score_min(self, thr_score_min: int) -> StringGroupSelection:
        """Create a new grouped selection with a different minimum score.

        Cached mapping-related grouped tables are reused. Edge and metric
        caches are not reused because they depend on the score threshold.

        Args:
            thr_score_min: Minimum `combined_score` required for retained
                STRING edges.

        Returns:
            A new grouped selection sharing cached grouped mapping state with
            the current selection.
        """
        return StringGroupSelection(
            dataset=self.dataset,
            _df_groups=self._df_groups,
            _df_input_ids=self._df_input_ids,
            thr_score_min=int(thr_score_min),
            _df_string_mapping=self._df_string_mapping,
            _df_unmapped_input_ids=self._df_unmapped_input_ids,
            _df_string_ids=self._df_string_ids,
        )

    def extract_string_mapping(self) -> pl.DataFrame:
        """Extract the grouped input-to-STRING mapping table.

        Returns:
            A materialized grouped mapping table with columns `GroupId`,
            `InputId`, `StringId`, and `MapSource`.

        Raises:
            ValueError: If the aliases file is missing required columns for the
                configured STRING version.
        """
        if self._df_string_mapping is None:
            alias_schema = self.dataset.alias_schema
            self._df_string_mapping = _extract_group_string_mapping_frame(
                file_aliases=self.dataset.file_aliases,
                df_input_ids=self._df_input_ids,
                source_rank_map=self.dataset.source_rank_map,
                col_string_id_aliases=alias_schema.col_string_id,
                has_source_aliases=alias_schema.has_source,
                version=self.dataset.version,
            )
        return self._df_string_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract grouped normalized input IDs that were not mapped.

        Returns:
            A materialized grouped table with columns `GroupId` and `InputId`
            for inputs that did not resolve through the aliases table.
        """
        if self._df_unmapped_input_ids is None:
            df_mapped_input_ids = self.extract_string_mapping().select(
                ["GroupId", "InputId"]
            )
            self._df_unmapped_input_ids = (
                self._df_input_ids.join(
                    df_mapped_input_ids.unique(),
                    on=["GroupId", "InputId"],
                    how="anti",
                )
                .select(["GroupId", "InputId"])
                .sort(["GroupId", "InputId"])
            )
        return self._df_unmapped_input_ids

    def extract_edges(self) -> pl.DataFrame:
        """Extract grouped STRING subnetworks induced by mapped STRING IDs.

        Returns:
            A materialized grouped edge table with columns `GroupId`,
            `StringIdA`, `StringIdB`, and `Score`.

        Raises:
            ValueError: If the links file is missing required columns for the
                configured STRING version.
        """
        if self._df_edges is None:
            self._df_edges = _extract_group_edges_frame(
                file_links=self.dataset.file_links,
                df_string_ids=self._extract_string_ids(),
                thr_score_min=self.thr_score_min,
                version=self.dataset.version,
            )
        return self._df_edges

    def extract_metrics(self) -> pl.DataFrame:
        """Extract grouped summary metrics derived from grouped outputs.

        Returns:
            A materialized grouped metrics table with one row per `GroupId` and
            count columns for inputs, mapped IDs, unmapped IDs, STRING IDs, and
            retained edges.
        """
        if self._df_metrics is None:
            self._df_metrics = _extract_group_metrics_frame(
                df_groups=self._df_groups,
                df_input_ids=self._df_input_ids,
                df_string_mapping=self.extract_string_mapping(),
                df_unmapped_input_ids=self.extract_unmapped_input_ids(),
                df_string_ids=self._extract_string_ids(),
                df_edges=self.extract_edges(),
            )
        return self._df_metrics

    def _extract_string_ids(self) -> pl.DataFrame:
        if self._df_string_ids is None:
            self._df_string_ids = (
                self.extract_string_mapping()
                .select(["GroupId", "StringId"])
                .unique(subset=["GroupId", "StringId"])
                .sort(["GroupId", "StringId"])
            )
        return self._df_string_ids
