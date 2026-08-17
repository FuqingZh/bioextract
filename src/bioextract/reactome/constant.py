from __future__ import annotations

from dataclasses import dataclass
from typing import TypedDict

import polars as pl
from polars._typing import SchemaDict

SCHEMA_VERSION = "reactome-mapping-v0.5"
SOURCE_SCHEMA_PROFILE = "reactome-mapping-files-v5"
MEDIA_TYPE_TSV = "text/tab-separated-values"
MEDIA_TYPE_ZIP = "application/zip"
ENTITY_COLUMN_MAPPING_REASON = "reactome_semantic_identifier_name"

MAPPING_LOWEST_LEVEL_ROLE = "uniprot_pathway_lowest_level"
MAPPING_ALL_LEVEL_ROLE = "uniprot_pathway_all_level"
PATHWAY_ROLE = "pathway"
RELATION_ROLE = "pathway_relation"
COMPLEX_PATHWAY_ROLE = "complex_pathway"
EWAS_PATHWAY_ROLE = "ewas_pathway"
PATHWAY_GENE_SET_ROLE = "pathway_gene_set"

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


@dataclass(frozen=True, slots=True)
class MappingRoleSpec:
    """Private contract for one official six-column mapping-family file."""

    role: str
    argument_name: str
    namespace: str
    target: str
    pathway_level: str | None
    filename: str
    source_column: str
    event_column: str
    name_column: str

    @property
    def raw_columns(self) -> tuple[str, ...]:
        return (
            self.source_column,
            self.event_column,
            "reactome_url",
            self.name_column,
            "evidence_code",
            "species",
        )

    @property
    def public_columns(self) -> tuple[str, ...]:
        return self.raw_columns


def _pathway_role(
    *,
    role: str,
    argument_name: str,
    namespace: str,
    filename: str,
    source_column: str,
) -> MappingRoleSpec:
    return MappingRoleSpec(
        role=role,
        argument_name=argument_name,
        namespace=namespace,
        target="pathway",
        pathway_level="lowest_level" if role.endswith("lowest_level") else "all_levels",
        filename=filename,
        source_column=source_column,
        event_column="reactome_pathway_id",
        name_column="pathway_name",
    )


def _reaction_role(
    *,
    role: str,
    argument_name: str,
    namespace: str,
    filename: str,
    source_column: str,
) -> MappingRoleSpec:
    return MappingRoleSpec(
        role=role,
        argument_name=argument_name,
        namespace=namespace,
        target="reaction",
        pathway_level=None,
        filename=filename,
        source_column=source_column,
        event_column="reactome_reaction_id",
        name_column="reaction_name",
    )


MAPPING_ROLE_SPECS: tuple[MappingRoleSpec, ...] = (
    _pathway_role(
        role=MAPPING_LOWEST_LEVEL_ROLE,
        argument_name="uniprot_mapping",
        namespace="uniprot",
        filename="UniProt2Reactome.txt",
        source_column="uniprot_id",
    ),
    _pathway_role(
        role=MAPPING_ALL_LEVEL_ROLE,
        argument_name="uniprot_all_levels",
        namespace="uniprot",
        filename="UniProt2Reactome_All_Levels.txt",
        source_column="uniprot_id",
    ),
    _reaction_role(
        role="uniprot_reaction",
        argument_name="uniprot_reactions",
        namespace="uniprot",
        filename="UniProt2ReactomeReactions.txt",
        source_column="uniprot_id",
    ),
    _pathway_role(
        role="ncbi_pathway_lowest_level",
        argument_name="ncbi_mapping",
        namespace="ncbi",
        filename="NCBI2Reactome.txt",
        source_column="ncbi_id",
    ),
    _pathway_role(
        role="ncbi_pathway_all_level",
        argument_name="ncbi_all_levels",
        namespace="ncbi",
        filename="NCBI2Reactome_All_Levels.txt",
        source_column="ncbi_id",
    ),
    _reaction_role(
        role="ncbi_reaction",
        argument_name="ncbi_reactions",
        namespace="ncbi",
        filename="NCBI2ReactomeReactions.txt",
        source_column="ncbi_id",
    ),
    _pathway_role(
        role="chebi_pathway_lowest_level",
        argument_name="chebi_mapping",
        namespace="chebi",
        filename="ChEBI2Reactome.txt",
        source_column="chebi_id",
    ),
    _pathway_role(
        role="chebi_pathway_all_level",
        argument_name="chebi_all_levels",
        namespace="chebi",
        filename="ChEBI2Reactome_All_Levels.txt",
        source_column="chebi_id",
    ),
    _reaction_role(
        role="chebi_reaction",
        argument_name="chebi_reactions",
        namespace="chebi",
        filename="ChEBI2ReactomeReactions.txt",
        source_column="chebi_id",
    ),
    _pathway_role(
        role="gtop_pathway_lowest_level",
        argument_name="gtop_mapping",
        namespace="gtop",
        filename="GtoP2Reactome.txt",
        source_column="gtop_id",
    ),
    _pathway_role(
        role="gtop_pathway_all_level",
        argument_name="gtop_all_levels",
        namespace="gtop",
        filename="GtoP2Reactome_All_Levels.txt",
        source_column="gtop_id",
    ),
    _reaction_role(
        role="gtop_reaction",
        argument_name="gtop_reactions",
        namespace="gtop",
        filename="GtoP2ReactomeReactions.txt",
        source_column="gtop_id",
    ),
)

MAPPING_ROLE_BY_ROLE = {spec.role: spec for spec in MAPPING_ROLE_SPECS}
MAPPING_ROLE_BY_ARGUMENT = {spec.argument_name: spec for spec in MAPPING_ROLE_SPECS}
MAPPING_ROLE_BY_DIMENSIONS = {
    (spec.namespace, spec.target, spec.pathway_level): spec
    for spec in MAPPING_ROLE_SPECS
}
MAPPING_OFFICIAL_FILENAMES = frozenset(spec.filename for spec in MAPPING_ROLE_SPECS)


class EntityRoleSpec(TypedDict):
    argument_name: str
    filename: str
    source_columns: tuple[str, str, str]
    public_columns: tuple[str, str, str]


class GmtSourceSpec(TypedDict):
    argument_name: str
    filename: str
    public_columns: tuple[str, str, str]


ENTITY_ROLE_SPECS: dict[str, EntityRoleSpec] = {
    COMPLEX_PATHWAY_ROLE: {
        "argument_name": "complex_pathways",
        "filename": "Complex_2_Pathway_human.txt",
        "source_columns": ("complex", "pathway", "top_level_pathway"),
        "public_columns": (
            "reactome_complex_id",
            "reactome_pathway_id",
            "top_level_reactome_pathway_id",
        ),
    },
    EWAS_PATHWAY_ROLE: {
        "argument_name": "ewas_pathways",
        "filename": "Ewas2Pathway_human.txt",
        "source_columns": ("ewas", "pathway", "top_level_pathway"),
        "public_columns": (
            "reactome_ewas_id",
            "reactome_pathway_id",
            "top_level_reactome_pathway_id",
        ),
    },
}
ENTITY_ROLE_BY_ARGUMENT = {
    value["argument_name"]: (role, value) for role, value in ENTITY_ROLE_SPECS.items()
}
ENTITY_ROLE_BY_ROLE = dict(ENTITY_ROLE_SPECS)

GMT_SOURCE_SPEC: GmtSourceSpec = {
    "argument_name": "pathway_gene_sets",
    "filename": "ReactomePathways.gmt.zip",
    "public_columns": (
        "reactome_pathway_id",
        "gene_set_name",
        "gene_symbol",
    ),
}

ASSET_SPECS = tuple(
    [(f"{spec.role}.parquet", "canonical", spec.role) for spec in MAPPING_ROLE_SPECS]
    + [
        ("pathway.parquet", "canonical", PATHWAY_ROLE),
        ("pathway_relation.parquet", "canonical", RELATION_ROLE),
        ("complex_pathway.parquet", "canonical", COMPLEX_PATHWAY_ROLE),
        ("ewas_pathway.parquet", "canonical", EWAS_PATHWAY_ROLE),
        ("pathway_gene_set.parquet", "canonical", PATHWAY_GENE_SET_ROLE),
    ]
)
