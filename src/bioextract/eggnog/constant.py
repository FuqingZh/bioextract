from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

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

SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}
