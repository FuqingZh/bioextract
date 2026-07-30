from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from pathlib import Path

import polars as pl

from .constant import (
    COLS_MAPPING,
    KIND_INPUT_ID_VALUES,
    SCHEMA_INTERPRO_ENTRY,
    SCHEMA_INTERPRO_MEMBER,
    SCHEMA_MAPPING,
    InterProInputIdKind,
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
            lf_mapping.drop("InterProType")
            .join(df_interpro_entry.lazy(), on="InterProId", how="left")
            .select(COLS_MAPPING)
        )
    if df_interpro_member.height > 0:
        lf_mapping = (
            lf_mapping.drop("MemberDb")
            .join(
                df_interpro_member.lazy(),
                on=["InterProId", "MemberDbId"],
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
                "UniProtId",
                "InterProId",
                "InterProName",
                "MemberDbId",
                "Start",
                "End",
            ],
            schema_overrides={
                "UniProtId": pl.String,
                "InterProId": pl.String,
                "InterProName": pl.String,
                "MemberDbId": pl.String,
                "Start": pl.Int64,
                "End": pl.Int64,
            },
            infer_schema=False,
            quote_char=None,
        )
        .with_columns(
            pl.lit(None, dtype=pl.String).alias("InterProType"),
            pl.lit(None, dtype=pl.String).alias("MemberDb"),
        )
        .select(COLS_MAPPING)
    )


def select_mapping_frame(
    file_protein2ipr: Path,
    df_input_ids: pl.DataFrame,
    *,
    kind_input_id: InterProInputIdKind,
    cols_group_id: tuple[str, ...],
    df_interpro_entry: pl.DataFrame,
    df_interpro_member: pl.DataFrame,
) -> pl.DataFrame:
    validate_kind_input_id(kind_input_id)
    if df_input_ids.height == 0:
        cols_out = list(cols_group_id) + ["InputId", "KindInputId"] + COLS_MAPPING
        return pl.DataFrame(
            schema={
                **dict.fromkeys(
                    list(cols_group_id) + ["InputId", "KindInputId"], pl.String
                ),
                **SCHEMA_MAPPING,
            }
        ).select(cols_out)

    input_ids = set(df_input_ids.get_column("InputId").to_list())
    return extract_mapping_frame(
        build_mapping_frame(
            iter_protein2ipr_records(file_protein2ipr, input_ids=input_ids),
            df_interpro_entry=df_interpro_entry,
            df_interpro_member=df_interpro_member,
        ),
        df_input_ids,
        kind_input_id=kind_input_id,
        cols_group_id=cols_group_id,
    )


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
            df_mapping.drop("InterProType")
            .join(df_interpro_entry, on="InterProId", how="left")
            .select(COLS_MAPPING)
        )
    if df_interpro_member.height > 0:
        df_mapping = (
            df_mapping.drop("MemberDb")
            .join(df_interpro_member, on=["InterProId", "MemberDbId"], how="left")
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
                "UniProtId": uniprot_id,
                "InterProId": interpro_id,
                "InterProName": interpro_name,
                "InterProType": None,
                "MemberDb": None,
                "MemberDbId": member_db_id,
                "Start": int(start),
                "End": int(end),
            }


def read_interpro_xml_frames(file_interpro_xml: Path | None) -> dict[str, pl.DataFrame]:
    if file_interpro_xml is None:
        return {
            "entry": pl.DataFrame(schema=SCHEMA_INTERPRO_ENTRY),
            "member": pl.DataFrame(schema=SCHEMA_INTERPRO_MEMBER),
        }

    entries: list[dict[str, str | None]] = []
    members: list[dict[str, str]] = []
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
                    "InterProId": interpro_id,
                    "InterProType": elem.attrib.get("type"),
                }
            )
            member_list = elem.find("member_list")
            if member_list is not None:
                for db_xref in member_list.findall("db_xref"):
                    member_db = db_xref.attrib.get("db")
                    member_db_id = db_xref.attrib.get("dbkey")
                    if member_db and member_db_id:
                        members.append(
                            {
                                "InterProId": interpro_id,
                                "MemberDbId": member_db_id,
                                "MemberDb": member_db,
                            }
                        )
            elem.clear()

    return {
        "entry": pl.DataFrame(entries, schema=SCHEMA_INTERPRO_ENTRY)
        .unique()
        .sort("InterProId")
        if entries
        else pl.DataFrame(schema=SCHEMA_INTERPRO_ENTRY),
        "member": pl.DataFrame(members, schema=SCHEMA_INTERPRO_MEMBER)
        .unique()
        .sort("InterProId", "MemberDbId")
        if members
        else pl.DataFrame(schema=SCHEMA_INTERPRO_MEMBER),
    }


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    kind_input_id: InterProInputIdKind,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    validate_kind_input_id(kind_input_id)
    cols_group = list(cols_group_id)
    cols_out = cols_group + ["InputId", "KindInputId"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="InputId",
            right_on="UniProtId",
            how="inner",
        )
        .with_columns(
            pl.col("InputId").alias("UniProtId"),
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
