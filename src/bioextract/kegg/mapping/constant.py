from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "kegg-mapping-v0.1"
MEDIA_TYPE_TSV = "text/tab-separated-values"

KIND_INPUT_ID_VALUES = ("uniprot", "ncbi_geneid", "kegg_gene")
KeggInputIdKind = Literal["uniprot", "ncbi_geneid", "kegg_gene"]

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

SCHEMA_UNMAPPED: SchemaDict = {"InputId": pl.String}
SCHEMA_GROUPS: SchemaDict = {"GroupId": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
}

ASSET_SPECS = (("mapping.parquet", "canonical", "mapping"),)
