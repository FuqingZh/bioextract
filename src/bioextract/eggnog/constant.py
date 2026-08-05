from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

NAMESPACE_VALUES = ("eggnog_protein",)
EggnogNamespace = Literal["eggnog_protein"]

COLS_MAPPING = [
    "name",
    "og",
    "level",
    "description",
    "COG_categories",
    "cog_category",
    "cog_class",
    "cog_name",
]

SCHEMA_MAPPING: SchemaDict = {
    "name": pl.String,
    "og": pl.String,
    "level": pl.String,
    "description": pl.String,
    "COG_categories": pl.String,
    "cog_category": pl.String,
    "cog_class": pl.String,
    "cog_name": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}
