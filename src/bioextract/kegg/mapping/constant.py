from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "kegg-mapping-v0.1"
MEDIA_TYPE_TSV = "text/tab-separated-values"

NAMESPACE_VALUES = ("uniprot", "ncbi_gene", "kegg_gene")
KEGGNamespace = Literal["uniprot", "ncbi_gene", "kegg_gene"]

COLS_MAPPING = [
    "OrganismCode",
    "KeggGeneId",
    "UniProtId",
    "NcbiGeneId",
    "KoId",
    "KeggPathwayId",
    "PathwayMapId",
    "GeneSymbol",
    "GeneDescription",
]

SCHEMA_MAPPING: SchemaDict = {
    "OrganismCode": pl.String,
    "KeggGeneId": pl.String,
    "UniProtId": pl.String,
    "NcbiGeneId": pl.String,
    "KoId": pl.String,
    "KeggPathwayId": pl.String,
    "PathwayMapId": pl.String,
    "GeneSymbol": pl.String,
    "GeneDescription": pl.String,
}

SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

ASSET_SPECS = (("mapping.parquet", "canonical", "mapping"),)
