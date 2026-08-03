from __future__ import annotations

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "uniprot-idmapping-duckdb-v1"
IDMAPPING_SOURCE_SCHEMA_PROFILE = "uniprot-idmapping-selected-22-column-v1"
SCHEMA_VERSION_EGGNOG_XREF = "uniprot-eggnog-xref-v0.1"
SCHEMA_VERSION_SUBCELLULAR_LOCATION = "uniprot-subcellular-location-v0.1"
MEDIA_TYPE_TSV_GZIP = "text/tab-separated-values+gzip"
MEDIA_TYPE_TSV = "text/tab-separated-values"
MEDIA_TYPE_PARQUET = "application/vnd.apache.parquet"
MEDIA_TYPE_PARQUET_DATASET = "application/vnd.apache.parquet.dataset"
MEDIA_TYPE_FLAT_FILE_GZIP = "application/gzip+uniprot-flat-file"
MEDIA_TYPE_FLAT_FILE = "application/vnd.uniprot.flat-file"

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

SCHEMA_MAPPING: SchemaDict = dict.fromkeys(COLS_IDMAPPING_SELECTED, pl.String)

REQUIRED_COLS_MAPPING = tuple(COLS_IDMAPPING_SELECTED)

COLS_EGGNOG_XREF = [
    "UniProtId",
    "PrimaryUniProtId",
    "IsPrimaryAccession",
    "EggnogOgId",
    "EggnogLevel",
    "SourceDb",
]

SCHEMA_EGGNOG_XREF: SchemaDict = {
    "UniProtId": pl.String,
    "PrimaryUniProtId": pl.String,
    "IsPrimaryAccession": pl.Boolean,
    "EggnogOgId": pl.String,
    "EggnogLevel": pl.String,
    "SourceDb": pl.String,
}

COLS_SUBCELLULAR_LOCATION = [
    "UniProtId",
    "PrimaryUniProtId",
    "UniProtEntryName",
    "GeneName",
    "ProteinName",
    "SubcellularLocation",
    "SubcellularLocationNote",
    "EvidenceCode",
    "EvidenceSource",
    "EvidenceId",
    "SourceDb",
]

SCHEMA_SUBCELLULAR_LOCATION: SchemaDict = {
    "UniProtId": pl.String,
    "PrimaryUniProtId": pl.String,
    "UniProtEntryName": pl.String,
    "GeneName": pl.String,
    "ProteinName": pl.String,
    "SubcellularLocation": pl.String,
    "SubcellularLocationNote": pl.String,
    "EvidenceCode": pl.String,
    "EvidenceSource": pl.String,
    "EvidenceId": pl.String,
    "SourceDb": pl.String,
}
