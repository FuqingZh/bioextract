from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "reactome-mapping-v0.1"
MEDIA_TYPE_TSV = "text/tab-separated-values"

COLS_MAPPING_RAW = [
    "UniProtId",
    "ReactomePathwayId",
    "ReactomeUrl",
    "PathwayName",
    "EvidenceCode",
    "Species",
]
COLS_PATHWAY_RAW = ["ReactomePathwayId", "PathwayName", "Species"]
COLS_RELATION_RAW = ["ParentReactomePathwayId", "ChildReactomePathwayId"]

SCHEMA_MAPPING_RAW: SchemaDict = {col: pl.String for col in COLS_MAPPING_RAW}
SCHEMA_PATHWAY_RAW: SchemaDict = {col: pl.String for col in COLS_PATHWAY_RAW}
SCHEMA_RELATION_RAW: SchemaDict = {col: pl.String for col in COLS_RELATION_RAW}

SCHEMA_GROUPS: SchemaDict = {"GroupId": pl.String}
SCHEMA_UNMAPPED: SchemaDict = {"InputId": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
}

ASSET_SPECS = (
    ("mapping.parquet", "canonical", "mapping"),
    ("pathway.parquet", "canonical", "pathway"),
    ("relation.parquet", "canonical", "relation"),
    ("term2gene.parquet", "derived", "term2gene"),
    ("term2name.parquet", "derived", "term2name"),
)
