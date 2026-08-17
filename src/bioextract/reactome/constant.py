from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "reactome-mapping-v0.2"
SOURCE_SCHEMA_PROFILE = "reactome-mapping-files-v2"
MEDIA_TYPE_TSV = "text/tab-separated-values"

MAPPING_LOWEST_LEVEL_ROLE = "uniprot_pathway_lowest_level"
MAPPING_ALL_LEVEL_ROLE = "uniprot_pathway_all_level"

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
    ("uniprot_pathway_lowest_level.parquet", "canonical", MAPPING_LOWEST_LEVEL_ROLE),
    ("uniprot_pathway_all_level.parquet", "canonical", MAPPING_ALL_LEVEL_ROLE),
    ("pathway.parquet", "canonical", "pathway"),
    ("pathway_relation.parquet", "canonical", "pathway_relation"),
)
