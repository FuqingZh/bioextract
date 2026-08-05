from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "wikipathways-gmt-v0.1"
MEDIA_TYPE_GMT = "application/gmt"

SCHEMA_PATHWAY: SchemaDict = {
    "wiki_pathways_id": pl.String,
    "pathway_name": pl.String,
    "species": pl.String,
    "collection": pl.String,
    "version": pl.String,
    "url": pl.String,
    "gene_count": pl.Int64,
}
SCHEMA_TERM2GENE: SchemaDict = {
    "wiki_pathways_id": pl.String,
    "gene_id": pl.String,
}
SCHEMA_TERM2NAME: SchemaDict = {
    "wiki_pathways_id": pl.String,
    "pathway_name": pl.String,
    "species": pl.String,
    "collection": pl.String,
    "version": pl.String,
    "url": pl.String,
}
SCHEMA_MAPPING: SchemaDict = {
    "input_id": pl.String,
    "gene_id": pl.String,
    "wiki_pathways_id": pl.String,
    "pathway_name": pl.String,
    "species": pl.String,
    "url": pl.String,
}
SCHEMA_GROUP_MAPPING: SchemaDict = {
    "group_id": pl.String,
    **SCHEMA_MAPPING,
}
SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

ASSET_SPECS = (
    ("pathway.parquet", "canonical", "pathway"),
    ("term2gene.parquet", "derived", "term2gene"),
    ("term2name.parquet", "derived", "term2name"),
)
