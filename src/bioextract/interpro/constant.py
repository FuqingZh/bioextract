from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "interpro-v1"
MEDIA_TYPE_XML_GZIP = "application/gzip+xml"
MEDIA_TYPE_TSV_GZIP = "application/gzip+tab-separated-values"

NAMESPACE_VALUES = ("uniprot",)
InterProNamespace = Literal["uniprot"]

COLS_MAPPING = [
    "uniprot_id",
    "interpro_id",
    "interpro_name",
    "interpro_type",
    "member_db",
    "member_db_id",
    "start",
    "end",
]

SCHEMA_MAPPING: SchemaDict = {
    "uniprot_id": pl.String,
    "interpro_id": pl.String,
    "interpro_name": pl.String,
    "interpro_type": pl.String,
    "member_db": pl.String,
    "member_db_id": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
}

SCHEMA_INTERPRO_ENTRY: SchemaDict = {
    "interpro_id": pl.String,
    "interpro_type": pl.String,
}

SCHEMA_INTERPRO_MEMBER: SchemaDict = {
    "interpro_id": pl.String,
    "member_db_id": pl.String,
    "member_db": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

ASSET_SPECS = (("mapping", "canonical", "mapping"),)
