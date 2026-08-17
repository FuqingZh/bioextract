from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from pathlib import Path

import polars as pl

from .constant import (
    COLS_MAPPING,
    NAMESPACE_VALUES,
    SCHEMA_INTERPRO_ENTRY,
    SCHEMA_INTERPRO_MEMBER,
    SCHEMA_MAPPING,
    InterProNamespace,
)


def read_mapping_frame(
    file_protein2ipr: Path,
    *,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> pl.DataFrame:
    return build_mapping_frame(
        iter_protein2ipr_records(file_protein2ipr),
        df_interpro_entry=df_interpro_entry,
        df_interpro_member=df_interpro_member,
    )


def scan_mapping_frame(
    file_protein2ipr: Path,
    *,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> pl.LazyFrame:
    lf_mapping = scan_protein2ipr_frame(file_protein2ipr)
    if df_interpro_entry.height > 0:
        lf_mapping = (
            lf_mapping.drop("interpro_type")
            .join(df_interpro_entry.lazy(), on="interpro_id", how="left")
            .select(COLS_MAPPING)
        )
    if df_interpro_member.height > 0:
        lf_mapping = (
            lf_mapping.drop("member_db")
            .join(
                df_interpro_member.lazy(),
                on=["interpro_id", "member_db_id"],
                how="left",
            )
            .select(COLS_MAPPING)
        )
    return lf_mapping.select(COLS_MAPPING)


def scan_protein2ipr_frame(file_protein2ipr: Path) -> pl.LazyFrame:
    return (
        pl.scan_csv(
            file_protein2ipr,
            separator="\t",
            has_header=False,
            new_columns=[
                "uniprot_id",
                "interpro_id",
                "interpro_name",
                "member_db_id",
                "start",
                "end",
            ],
            schema_overrides={
                "uniprot_id": pl.String,
                "interpro_id": pl.String,
                "interpro_name": pl.String,
                "member_db_id": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
            },
            infer_schema=False,
            quote_char=None,
        )
        .with_columns(
            pl.lit(None, dtype=pl.String).alias("interpro_type"),
            pl.lit(None, dtype=pl.String).alias("member_db"),
        )
        .select(COLS_MAPPING)
    )


def validate_mapping_xml_relationships(
    file_protein2ipr: Path,
    *,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> None:
    """Require one complete XML explanation for every mapping relationship."""
    mapping_keys = (
        scan_protein2ipr_frame(file_protein2ipr)
        .select("interpro_id", "member_db_id")
        .unique()
        .collect(engine="streaming")
    )
    entry_values: dict[str, set[str]] = {}
    entry_row_counts: dict[str, int] = {}
    invalid_entry_keys: set[str] = set()
    entry_keys: set[str] = set()
    for interpro_id, interpro_type in df_interpro_entry.iter_rows():
        entry_keys.add(interpro_id)
        entry_row_counts[interpro_id] = entry_row_counts.get(interpro_id, 0) + 1
        if interpro_type is not None and interpro_type.strip():
            entry_values.setdefault(interpro_id, set()).add(interpro_type.strip())
        else:
            invalid_entry_keys.add(interpro_id)
    member_values: dict[tuple[str, str], set[str]] = {}
    member_row_counts: dict[tuple[str, str], int] = {}
    invalid_member_keys: set[tuple[str, str]] = set()
    member_keys: set[tuple[str, str]] = set()
    for interpro_id, member_db_id, member_db in df_interpro_member.iter_rows():
        key = (interpro_id, member_db_id)
        member_keys.add(key)
        member_row_counts[key] = member_row_counts.get(key, 0) + 1
        if member_db is not None and member_db.strip():
            member_values.setdefault(key, set()).add(member_db.strip())
        else:
            invalid_member_keys.add(key)

    for interpro_id, member_db_id in mapping_keys.iter_rows():
        if interpro_id not in entry_keys:
            raise ValueError(
                "InterPro mapping relationship is absent from XML entry metadata: "
                f"{interpro_id}"
            )
        types = entry_values.get(interpro_id, set())
        if (
            entry_row_counts[interpro_id] != 1
            or interpro_id in invalid_entry_keys
            or len(types) != 1
        ):
            raise ValueError(
                "InterPro mapping entry is not uniquely explained by one non-empty "
                f"XML InterProType: {interpro_id}"
            )
        member_key = (interpro_id, member_db_id)
        if member_key not in member_keys:
            raise ValueError(
                "InterPro mapping relationship is absent from XML member metadata: "
                f"{interpro_id}/{member_db_id}"
            )
        databases = member_values.get(member_key, set())
        if (
            member_row_counts[member_key] != 1
            or member_key in invalid_member_keys
            or len(databases) != 1
        ):
            raise ValueError(
                "InterPro mapping member is not uniquely explained by one non-empty "
                f"XML MemberDb: {interpro_id}/{member_db_id}"
            )


def select_mapping_frame(
    file_protein2ipr: Path,
    df_input_ids: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
    namespace: InterProNamespace,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> pl.DataFrame:
    validate_namespace(namespace)
    if df_input_ids.height == 0:
        cols_group = ["group_id"] if df_group_membership is not None else []
        cols_out = cols_group + ["input_id", "input_namespace"] + COLS_MAPPING
        return pl.DataFrame(
            schema={
                **dict.fromkeys(
                    cols_group + ["input_id", "input_namespace"], pl.String
                ),
                **SCHEMA_MAPPING,
            }
        ).select(cols_out)

    input_ids = set(df_input_ids.get_column("input_id").to_list())
    df_selected = extract_mapping_frame(
        build_mapping_frame(
            iter_protein2ipr_records(file_protein2ipr, input_ids=input_ids),
            df_interpro_entry=df_interpro_entry,
            df_interpro_member=df_interpro_member,
        ),
        df_input_ids,
        namespace=namespace,
    )
    return _expand_group_membership(df_selected, df_group_membership)


def build_mapping_frame(
    records: Iterable[dict[str, str | int | None]],
    *,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> pl.DataFrame:
    rows = list(records)
    if not rows:
        return pl.DataFrame(schema=SCHEMA_MAPPING)

    df_mapping = pl.DataFrame(rows, schema=SCHEMA_MAPPING)
    if df_interpro_entry.height > 0:
        df_mapping = (
            df_mapping.drop("interpro_type")
            .join(df_interpro_entry, on="interpro_id", how="left")
            .select(COLS_MAPPING)
        )
    if df_interpro_member.height > 0:
        df_mapping = (
            df_mapping.drop("member_db")
            .join(df_interpro_member, on=["interpro_id", "member_db_id"], how="left")
            .select(COLS_MAPPING)
        )
    return df_mapping.unique().sort(COLS_MAPPING)


def iter_protein2ipr_records(
    file_protein2ipr: Path,
    *,
    input_ids: set[str] | None = None,
) -> Iterable[dict[str, str | int | None]]:
    with gzip.open(
        file_protein2ipr, "rt", encoding="utf-8", errors="replace"
    ) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 6:
                raise ValueError(
                    "InterPro protein2ipr line must contain 6 columns: "
                    f"path={file_protein2ipr}, line={line_number}, value={line!r}"
                )
            uniprot_id, interpro_id, interpro_name, member_db_id, start, end = fields
            if input_ids is not None and uniprot_id not in input_ids:
                continue
            yield {
                "uniprot_id": uniprot_id,
                "interpro_id": interpro_id,
                "interpro_name": interpro_name,
                "interpro_type": None,
                "member_db": None,
                "member_db_id": member_db_id,
                "start": int(start),
                "end": int(end),
            }


def iter_mapping_frames(
    file_protein2ipr: Path,
    *,
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
    input_ids: set[str] | None = None,
    _batch_size: int = 10_000,
) -> Iterable[pl.DataFrame]:
    """Yield bounded mapping frames while parsing a protein2ipr source."""
    rows: list[dict[str, str | int | None]] = []
    for record in iter_protein2ipr_records(file_protein2ipr, input_ids=input_ids):
        rows.append(record)
        if len(rows) >= _batch_size:
            yield build_mapping_frame(
                rows,
                df_interpro_entry=df_interpro_entry,
                df_interpro_member=df_interpro_member,
            )
            rows = []
    if rows:
        yield build_mapping_frame(
            rows,
            df_interpro_entry=df_interpro_entry,
            df_interpro_member=df_interpro_member,
        )


def read_interpro_xml_frames(file_interpro_xml: Path | None) -> dict[str, pl.DataFrame]:
    if file_interpro_xml is None:
        return {
            "entry": pl.DataFrame(schema=SCHEMA_INTERPRO_ENTRY),
            "member": pl.DataFrame(schema=SCHEMA_INTERPRO_MEMBER),
        }

    entries: list[dict[str, str | None]] = []
    members: list[dict[str, str | None]] = []
    with gzip.open(file_interpro_xml, "rb") as handle:
        for _event, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag != "interpro":
                continue
            interpro_id = elem.attrib.get("id")
            if interpro_id is None:
                elem.clear()
                continue
            entries.append(
                {
                    "interpro_id": interpro_id,
                    "interpro_type": elem.attrib.get("type"),
                }
            )
            member_list = elem.find("member_list")
            if member_list is not None:
                for db_xref in member_list.findall("db_xref"):
                    member_db = db_xref.attrib.get("db")
                    member_db_id = db_xref.attrib.get("dbkey")
                    if member_db_id:
                        members.append(
                            {
                                "interpro_id": interpro_id,
                                "member_db_id": member_db_id,
                                "member_db": member_db,
                            }
                        )
            elem.clear()

    return {
        "entry": pl.DataFrame(entries, schema=SCHEMA_INTERPRO_ENTRY)
        .unique()
        .sort("interpro_id")
        if entries
        else pl.DataFrame(schema=SCHEMA_INTERPRO_ENTRY),
        "member": pl.DataFrame(members, schema=SCHEMA_INTERPRO_MEMBER)
        .unique()
        .sort("interpro_id", "member_db_id")
        if members
        else pl.DataFrame(schema=SCHEMA_INTERPRO_MEMBER),
    }


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    namespace: InterProNamespace,
) -> pl.DataFrame:
    validate_namespace(namespace)
    cols_out = ["input_id", "input_namespace"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="input_id",
            right_on="uniprot_id",
            how="inner",
        )
        .with_columns(
            pl.col("input_id").alias("uniprot_id"),
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
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    df_mapped_input_ids = df_mapping.select("input_id").unique().sort("input_id")
    df_unmatched = (
        df_input_ids.join(df_mapped_input_ids, on="input_id", how="anti")
        .select("input_id")
        .sort("input_id")
    )
    if df_group_membership is None:
        return df_unmatched
    return (
        df_group_membership.join(df_unmatched, on="input_id", how="inner")
        .select("group_id", "input_id")
        .sort("group_id", "input_id")
    )


def _expand_group_membership(
    df_mapping: pl.DataFrame,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    if df_group_membership is None:
        return df_mapping
    cols_out = ["group_id", *df_mapping.columns]
    return (
        df_group_membership.join(df_mapping, on="input_id", how="inner")
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def validate_namespace(namespace: str) -> None:
    if namespace not in NAMESPACE_VALUES:
        raise ValueError(
            "namespace must be one of: "
            f"{', '.join(NAMESPACE_VALUES)}; got {namespace!r}"
        )
