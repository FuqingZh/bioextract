from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "kegg-mapping-v0.1"
MEDIA_TYPE_TSV = "text/tab-separated-values"

NAMESPACE_VALUES = ("uniprot", "ncbi_gene", "kegg_gene")
KEGGNamespace = Literal["uniprot", "ncbi_gene", "kegg_gene"]

COLS_MAPPING = [
    "organism_code",
    "kegg_gene_id",
    "uniprot_id",
    "ncbi_gene_id",
    "ko_id",
    "kegg_pathway_id",
    "pathway_map_id",
    "gene_symbol",
    "gene_description",
]

SCHEMA_MAPPING: SchemaDict = {
    "organism_code": pl.String,
    "kegg_gene_id": pl.String,
    "uniprot_id": pl.String,
    "ncbi_gene_id": pl.String,
    "ko_id": pl.String,
    "kegg_pathway_id": pl.String,
    "pathway_map_id": pl.String,
    "gene_symbol": pl.String,
    "gene_description": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

ASSET_SPECS = (("mapping.parquet", "canonical", "mapping"),)
