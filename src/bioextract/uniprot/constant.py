from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "uniprot-idmapping-selected-v0.1"
MEDIA_TYPE_TSV_GZIP = "text/tab-separated-values+gzip"
MEDIA_TYPE_TSV = "text/tab-separated-values"
MEDIA_TYPE_PARQUET = "application/vnd.apache.parquet"
MEDIA_TYPE_PARQUET_DATASET = "application/vnd.apache.parquet.dataset"

COLS_IDMAPPING_SELECTED = [
    "UniProtId",
    "UniProtEntryName",
    "GeneId",
    "RefSeq",
    "GI",
    "PDB",
    "GO",
    "UniRef100",
    "UniRef90",
    "UniRef50",
    "UniParc",
    "PIR",
    "TaxId",
    "MIM",
    "UniGene",
    "PubMed",
    "EMBL",
    "EMBLCDS",
    "Ensembl",
    "EnsemblTranscript",
    "EnsemblProtein",
    "AdditionalPubMed",
]

SCHEMA_MAPPING: SchemaDict = {col: pl.String for col in COLS_IDMAPPING_SELECTED}

REQUIRED_COLS_MAPPING = tuple(COLS_IDMAPPING_SELECTED)
