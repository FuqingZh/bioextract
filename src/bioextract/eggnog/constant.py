from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "eggnog-mapping-v0.1"
MEDIA_TYPE_SQLITE_GZIP = "application/gzip+sqlite3"
MEDIA_TYPE_SQLITE = "application/vnd.sqlite3"
MEDIA_TYPE_TSV = "text/tab-separated-values"

NAMESPACE_VALUES = ("eggnog_protein",)
EggnogNamespace = Literal["eggnog_protein"]

COLS_MAPPING = [
    "EggnogProteinId",
    "EggnogOgId",
    "EggnogLevel",
    "CogCategory",
    "CogClass",
    "CogName",
    "OgDescription",
]

SCHEMA_MAPPING: SchemaDict = {
    "EggnogProteinId": pl.String,
    "EggnogOgId": pl.String,
    "EggnogLevel": pl.String,
    "CogCategory": pl.String,
    "CogClass": pl.String,
    "CogName": pl.String,
    "OgDescription": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"InputId": pl.String}
SCHEMA_GROUPS: SchemaDict = {"GroupId": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
}

ASSET_SPECS = (("mapping.parquet", "canonical", "mapping"),)
