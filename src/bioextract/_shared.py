import csv
import re
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass
from typing import Protocol, TextIO

import polars as pl
from polars._typing import SchemaDict

_RE_UNIPROT_SELECTION = re.compile(r"^(?:sp|tr)\|([^|]+)\|([^|]+)$")


@dataclass(frozen=True, slots=True)
class GroupInputFrames:
    """Carry normalized groups, membership, and globally unique input IDs."""

    df_groups: pl.DataFrame
    df_group_membership: pl.DataFrame
    df_input_ids: pl.DataFrame


class RowWriter(Protocol):
    """Minimal row-writer contract shared by streaming TSV producers."""

    def writerow(self, row: Iterable[object], /) -> object: ...


def create_tsv_writer(handle: TextIO) -> RowWriter:
    """Create a tab-delimited writer with deterministic LF row endings."""
    return csv.writer(handle, delimiter="\t", lineterminator="\n")


def normalize_uniprot_selection_id(value: object) -> str:
    """Normalize one caller-supplied UniProt representation.

    Plain non-pipe text is preserved after trimming. Pipe-bearing text must be
    one complete sp|accession|entry_name or tr|accession|entry_name value;
    malformed or non-UniProt pipe forms are rejected instead of becoming
    unmatched identifiers.

    Examples:
        >>> normalize_uniprot_selection_id(" sp|P04637|P53_HUMAN ")
        'P04637'
        >>> normalize_uniprot_selection_id(" P04637 ")
        'P04637'
        >>> normalize_uniprot_selection_id("db|P04637|P53_HUMAN")
        Traceback (most recent call last):
        ...
        ValueError: Invalid UniProt pipe-form identifier: 'db|P04637|P53_HUMAN'
    """
    normalized = str(value).strip()
    if "|" not in normalized:
        return normalized
    match = _RE_UNIPROT_SELECTION.fullmatch(normalized)
    if match is None:
        raise ValueError(f"Invalid UniProt pipe-form identifier: {normalized!r}")
    accession = match.group(1).strip()
    entry_name = match.group(2).strip()
    if not accession or not entry_name:
        raise ValueError(f"Invalid UniProt pipe-form identifier: {normalized!r}")
    return accession


def validate_group_ids(group_ids: list[str]) -> None:
    group_ids_seen: set[str] = set()
    for group_id in group_ids:
        if not group_id:
            raise ValueError("group_id must be a non-empty string after normalization")
        if group_id in group_ids_seen:
            raise ValueError(
                f"group_id values must be unique after normalization: {group_id!r}"
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
        input_ids: Raw identifiers to trim.
        schema_unmapped: Output schema containing the `input_id` column.

    Returns:
        A table of non-empty, unique trimmed IDs sorted by `input_id`.
    """
    ids_normalized: list[str] = []
    for input_id in input_ids:
        if input_id_normalized := str(input_id).strip():
            ids_normalized.append(input_id_normalized)

    if not ids_normalized:
        return pl.DataFrame(schema=schema_unmapped)

    return (
        pl.DataFrame({"input_id": ids_normalized}, schema=schema_unmapped)
        .unique(subset=["input_id"])
        .sort("input_id")
    )


def create_group_input_frames(
    ids_by_group: Mapping[str, Iterable[str]],
    *,
    schema_groups: SchemaDict,
    schema_group_input_ids: SchemaDict,
) -> GroupInputFrames:
    """Build canonical group and grouped-input tables.

    Group labels and input IDs are stripped. Group labels must remain non-empty
    and unique. Membership rows are deduplicated within each group, while the
    input table contains one globally unique row per trimmed ID. Groups with no
    retained IDs remain present in the group registry.

    Args:
        ids_by_group: Mapping of raw group labels to raw identifiers.
        schema_groups: Output schema containing `group_id`.
        schema_group_input_ids: Output schema containing `group_id` and
            `input_id`.

    Returns:
        Sorted group-registry, group-membership, and unique-input tables.

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
            if input_id_normalized := str(input_id).strip():
                group_ids_col.append(group_id)
                input_ids_col.append(input_id_normalized)

    validate_group_ids(group_ids_normalized)

    if group_ids_normalized:
        df_groups = pl.DataFrame(
            {"group_id": group_ids_normalized},
            schema=schema_groups,
        ).sort("group_id")
    else:
        df_groups = pl.DataFrame(schema=schema_groups)

    if not input_ids_col:
        df_group_membership = pl.DataFrame(schema=schema_group_input_ids)
        return GroupInputFrames(
            df_groups=df_groups,
            df_group_membership=df_group_membership,
            df_input_ids=df_group_membership.select("input_id"),
        )

    df_group_membership = (
        pl.DataFrame(
            {"group_id": group_ids_col, "input_id": input_ids_col},
            schema=schema_group_input_ids,
        )
        .unique()
        .sort("group_id", "input_id")
    )
    df_input_ids = df_group_membership.select("input_id").unique().sort("input_id")
    return GroupInputFrames(
        df_groups=df_groups,
        df_group_membership=df_group_membership,
        df_input_ids=df_input_ids,
    )
