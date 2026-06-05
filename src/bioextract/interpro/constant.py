from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "interpro-mapping-v0.1"
MEDIA_TYPE_XML_GZIP = "application/gzip+xml"
MEDIA_TYPE_TSV_GZIP = "application/gzip+tab-separated-values"

KIND_INPUT_ID_VALUES = ("uniprot",)
InterProInputIdKind = Literal["uniprot"]

COLS_MAPPING = [
    "UniProtId",
    "InterProId",
    "InterProName",
    "InterProType",
    "MemberDb",
    "MemberDbId",
    "Start",
    "End",
]

SCHEMA_MAPPING: SchemaDict = {
    "UniProtId": pl.String,
    "InterProId": pl.String,
    "InterProName": pl.String,
    "InterProType": pl.String,
    "MemberDb": pl.String,
    "MemberDbId": pl.String,
    "Start": pl.Int64,
    "End": pl.Int64,
}

SCHEMA_INTERPRO_ENTRY: SchemaDict = {
    "InterProId": pl.String,
    "InterProType": pl.String,
}

SCHEMA_INTERPRO_MEMBER: SchemaDict = {
    "InterProId": pl.String,
    "MemberDbId": pl.String,
    "MemberDb": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"InputId": pl.String}
SCHEMA_GROUPS: SchemaDict = {"GroupId": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
}

ASSET_SPECS = (("mapping.parquet", "canonical", "mapping"),)
