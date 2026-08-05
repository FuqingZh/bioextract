from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "reactome-mapping-v0.1"
MEDIA_TYPE_TSV = "text/tab-separated-values"

COLS_MAPPING_RAW = [
    "uniprot_id",
    "reactome_pathway_id",
    "reactome_url",
    "pathway_name",
    "evidence_code",
    "species",
]
COLS_PATHWAY_RAW = ["reactome_pathway_id", "pathway_name", "species"]
COLS_RELATION_RAW = [
    "parent_reactome_pathway_id",
    "child_reactome_pathway_id",
]

SCHEMA_MAPPING_RAW: SchemaDict = dict.fromkeys(COLS_MAPPING_RAW, pl.String)
SCHEMA_PATHWAY_RAW: SchemaDict = dict.fromkeys(COLS_PATHWAY_RAW, pl.String)
SCHEMA_RELATION_RAW: SchemaDict = dict.fromkeys(COLS_RELATION_RAW, pl.String)

SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

ASSET_SPECS = (
    ("mapping.parquet", "canonical", "mapping"),
    ("pathway.parquet", "canonical", "pathway"),
    ("relation.parquet", "canonical", "relation"),
    ("term2gene.parquet", "derived", "term2gene"),
    ("term2name.parquet", "derived", "term2name"),
)
