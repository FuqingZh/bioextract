import csv
import re
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass
from typing import Protocol, TextIO

import polars as pl
from polars._typing import SchemaDict

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")


@dataclass(frozen=True, slots=True)
class GroupInputFrames:
    """Carry a normalized group registry and its deduplicated input-ID rows."""

    df_groups: pl.DataFrame
    df_input_ids: pl.DataFrame


class RowWriter(Protocol):
    """Minimal row-writer contract shared by streaming TSV producers."""

    def writerow(self, row: Iterable[object], /) -> object: ...


def create_tsv_writer(handle: TextIO) -> RowWriter:
    """Create a tab-delimited writer with deterministic LF row endings."""
    return csv.writer(handle, delimiter="\t", lineterminator="\n")


def normalize_input_id(value: str) -> str:
    """Normalize an identifier shared by all selection APIs.

    Surrounding whitespace is removed. UniProt pipe-style identifiers are
    reduced to the accession between the first two pipes; other identifiers
    are returned unchanged after trimming.

    Examples:
        >>> normalize_input_id(" sp|P04637|P53_HUMAN ")
        'P04637'
        >>> normalize_input_id(" TP53 ")
        'TP53'
    """
    value = value.strip()
    if (match_pipe := RE_UNIPROT_PIPE.match(value)) is not None:
        return match_pipe.group(1).strip()
    return value


def validate_group_ids(group_ids: list[str]) -> None:
    group_ids_seen: set[str] = set()
    for group_id in group_ids:
        if not group_id:
            raise ValueError("GroupId must be a non-empty string after normalization")
        if group_id in group_ids_seen:
            raise ValueError(
                f"GroupId values must be unique after normalization: {group_id!r}"
            )
        group_ids_seen.add(group_id)


def validate_required_cols(
    cols_available: Collection[str], cols_required: Collection[str], context: str
) -> None:
    cols_missing = set(cols_required) - set(cols_available)
    if cols_missing:
        raise ValueError(
            f"{context} is missing required columns: "
            f"{sorted(cols_missing)}; available={cols_available}"
        )


def create_input_id_frame(
    input_ids: Iterable[str],
    *,
    schema_unmapped: SchemaDict,
) -> pl.DataFrame:
    """Build the canonical single-selection input table.

    Args:
        input_ids: Raw identifiers to normalize.
        schema_unmapped: Output schema containing the `InputId` column.

    Returns:
        A table of non-empty, unique normalized IDs sorted by `InputId`.
    """
    ids_normalized: list[str] = []
    for input_id in input_ids:
        if input_id_normalized := normalize_input_id(str(input_id)):
            ids_normalized.append(input_id_normalized)

    if not ids_normalized:
        return pl.DataFrame(schema=schema_unmapped)

    return (
        pl.DataFrame({"InputId": ids_normalized}, schema=schema_unmapped)
        .unique(subset=["InputId"])
        .sort("InputId")
    )


def create_group_input_frames(
    ids_by_group: Mapping[str, Iterable[str]],
    *,
    schema_groups: SchemaDict,
    schema_group_input_ids: SchemaDict,
) -> GroupInputFrames:
    """Build canonical group and grouped-input tables.

    Group labels are stripped and must remain non-empty and unique. Input IDs
    use `normalize_input_id()` and are deduplicated within each group. Groups
    with no retained IDs remain present in the group registry.

    Args:
        ids_by_group: Mapping of raw group labels to raw identifiers.
        schema_groups: Output schema containing `GroupId`.
        schema_group_input_ids: Output schema containing `GroupId` and
            `InputId`.

    Returns:
        Sorted group-registry and grouped-input tables.

    Raises:
        ValueError: If a normalized group label is empty or duplicated.
    """
    group_ids_normalized: list[str] = []
    group_ids_col: list[str] = []
    input_ids_col: list[str] = []

    for group_id_raw, ids in ids_by_group.items():
        group_id = str(group_id_raw).strip()
        group_ids_normalized.append(group_id)
        for input_id in ids:
            if input_id_normalized := normalize_input_id(str(input_id)):
                group_ids_col.append(group_id)
                input_ids_col.append(input_id_normalized)

    validate_group_ids(group_ids_normalized)

    if group_ids_normalized:
        df_groups = pl.DataFrame(
            {"GroupId": group_ids_normalized},
            schema=schema_groups,
        ).sort("GroupId")
    else:
        df_groups = pl.DataFrame(schema=schema_groups)

    if not input_ids_col:
        return GroupInputFrames(
            df_groups=df_groups,
            df_input_ids=pl.DataFrame(schema=schema_group_input_ids),
        )

    df_group_input_ids = (
        pl.DataFrame(
            {"GroupId": group_ids_col, "InputId": input_ids_col},
            schema=schema_group_input_ids,
        )
        .unique()
        .sort("GroupId", "InputId")
    )
    return GroupInputFrames(df_groups=df_groups, df_input_ids=df_group_input_ids)
