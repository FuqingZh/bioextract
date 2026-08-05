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
    "uniprot_id",
    "uniprot_entry_name",
    "gene_id",
    "refseq",
    "gi",
    "pdb",
    "go",
    "uniref100",
    "uniref90",
    "uniref50",
    "uniparc",
    "pir",
    "tax_id",
    "mim",
    "unigene",
    "pubmed",
    "embl",
    "embl_cds",
    "ensembl",
    "ensembl_transcript",
    "ensembl_protein",
    "additional_pubmed",
]

SCHEMA_MAPPING: SchemaDict = dict.fromkeys(COLS_IDMAPPING_SELECTED, pl.String)

REQUIRED_COLS_MAPPING = tuple(COLS_IDMAPPING_SELECTED)

COLS_EGGNOG_XREF = [
    "uniprot_id",
    "primary_uniprot_id",
    "is_primary_accession",
    "eggnog_og_id",
    "eggnog_level",
    "source_db",
]

SCHEMA_EGGNOG_XREF: SchemaDict = {
    "uniprot_id": pl.String,
    "primary_uniprot_id": pl.String,
    "is_primary_accession": pl.Boolean,
    "eggnog_og_id": pl.String,
    "eggnog_level": pl.String,
    "source_db": pl.String,
}

COLS_SUBCELLULAR_LOCATION = [
    "uniprot_id",
    "primary_uniprot_id",
    "uniprot_entry_name",
    "gene_name",
    "protein_name",
    "subcellular_location",
    "subcellular_location_note",
    "evidence_code",
    "evidence_source",
    "evidence_id",
    "source_db",
]

SCHEMA_SUBCELLULAR_LOCATION: SchemaDict = {
    "uniprot_id": pl.String,
    "primary_uniprot_id": pl.String,
    "uniprot_entry_name": pl.String,
    "gene_name": pl.String,
    "protein_name": pl.String,
    "subcellular_location": pl.String,
    "subcellular_location_note": pl.String,
    "evidence_code": pl.String,
    "evidence_source": pl.String,
    "evidence_id": pl.String,
    "source_db": pl.String,
}
