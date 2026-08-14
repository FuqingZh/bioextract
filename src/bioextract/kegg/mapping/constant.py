from __future__ import annotations

from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "kegg-mapping-v1.0"
SOURCE_SCHEMA_PROFILE = "kegg-organism-mapping-files-v2"
MEDIA_TYPE_TSV = "text/tab-separated-values"

NAMESPACE_VALUES = ("uniprot", "ncbi_gene", "kegg_gene")
KEGGNamespace = Literal["uniprot", "ncbi_gene", "kegg_gene"]

ORGANISM_ROLE_FILENAMES = {
    "gene_list": "gene_list.tsv",
    "uniprot_conversion": "conv_uniprot.tsv",
    "ncbi_gene_conversion": "conv_ncbi_geneid.tsv",
    "gene_ko": "gene_ko.tsv",
    "gene_pathway": "gene_pathway.tsv",
}
CAPABILITY_NAMES = (
    *ORGANISM_ROLE_FILENAMES,
    "organism_list",
    "ko_pathway",
)

SCHEMA_ORGANISM: SchemaDict = {
    "organism_code": pl.String,
    "genome_id": pl.String,
    "organism_name": pl.String,
    "taxonomy_lineage": pl.List(pl.String),
}

SCHEMA_GENE_ANNOTATION: SchemaDict = {
    "organism_code": pl.String,
    "kegg_gene_id": pl.String,
    "gene_type": pl.String,
    "genomic_position": pl.String,
    "gene_symbol": pl.String,
    "gene_aliases": pl.List(pl.String),
    "gene_description": pl.String,
    "uniprot_mappings": pl.List(pl.Struct({"uniprot_id": pl.String})),
    "ncbi_gene_mappings": pl.List(pl.Struct({"ncbi_gene_id": pl.String})),
    "ko_mappings": pl.List(pl.Struct({"ko_id": pl.String})),
    "pathway_mappings": pl.List(
        pl.Struct(
            {
                "kegg_pathway_id": pl.String,
                "pathway_map_id": pl.String,
            }
        )
    ),
}

SCHEMA_KO_ANNOTATION: SchemaDict = {
    "ko_id": pl.String,
    "pathway_mappings": pl.List(
        pl.Struct(
            {
                "kegg_pathway_id": pl.String,
                "pathway_namespace": pl.String,
                "pathway_map_id": pl.String,
            }
        )
    ),
}

SCHEMA_GENE_PATHWAY: SchemaDict = {
    "organism_code": pl.String,
    "kegg_gene_id": pl.String,
    "pathway_mappings": SCHEMA_GENE_ANNOTATION["pathway_mappings"],
}
SCHEMA_GENE_PATHWAY_VIA_KO: SchemaDict = {
    "organism_code": pl.String,
    "kegg_gene_id": pl.String,
    "pathway_mappings": pl.List(
        pl.Struct(
            {
                "ko_id": pl.String,
                "kegg_pathway_id": pl.String,
                "pathway_namespace": pl.String,
                "pathway_map_id": pl.String,
            }
        )
    ),
}
SCHEMA_KO_PATHWAY: SchemaDict = dict(SCHEMA_KO_ANNOTATION)
SCHEMA_MATCH: SchemaDict = {
    "input_id": pl.String,
    "input_namespace": pl.String,
    "organism_code": pl.String,
    "kegg_gene_id": pl.String,
}
SCHEMA_UNMAPPED: SchemaDict = {"input_id": pl.String}
SCHEMA_GROUPS: SchemaDict = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
}

TABLE_SCHEMAS = {
    "organism": SCHEMA_ORGANISM,
    "gene_annotation": SCHEMA_GENE_ANNOTATION,
    "ko_annotation": SCHEMA_KO_ANNOTATION,
}
