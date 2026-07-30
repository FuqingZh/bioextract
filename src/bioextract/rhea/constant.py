"""Stable source and schema constants for Rhea extraction."""

from __future__ import annotations

from typing import Literal

SCHEMA_VERSION = "rhea-duckdb-v1"

RheaNamespace = Literal[
    "rhea",
    "chebi",
    "uniprot",
    "ec",
    "go",
    "ecocyc",
    "kegg_reaction",
    "macie",
    "metacyc",
    "reactome",
]

SOURCE_BASENAMES: dict[str, tuple[str, ...]] = {
    "license": ("LICENSE.txt",),
    "chebi_names": ("chebiId_name.tsv", "chebiId_name.tsv.gz"),
    "chebi_ph7_3_mapping": (
        "chebi_pH7_3_mapping.tsv",
        "chebi_pH7_3_mapping.tsv.gz",
    ),
    "directions": ("rhea-directions.tsv", "rhea-directions.tsv.gz"),
    "obsoletes": ("rhea-obsoletes.tsv", "rhea-obsoletes.tsv.gz"),
    "reaction_smiles": (
        "rhea-reaction-smiles.tsv",
        "rhea-reaction-smiles.tsv.gz",
    ),
    "relationships": (
        "rhea-relationships.tsv",
        "rhea-relationships.tsv.gz",
    ),
    "release_properties": ("rhea-release.properties",),
    "rdf": ("rhea.rdf", "rhea.rdf.gz"),
    "sdf": ("rhea.sdf", "rhea.sdf.gz"),
    "ec": ("rhea2ec.tsv", "rhea2ec.tsv.gz"),
    "go": ("rhea2go.tsv", "rhea2go.tsv.gz"),
    "uniprot_sprot": (
        "rhea2uniprot_sprot.tsv",
        "rhea2uniprot_sprot.tsv.gz",
    ),
    "uniprot_trembl": (
        "rhea2uniprot_trembl.tsv",
        "rhea2uniprot_trembl.tsv.gz",
    ),
    "xrefs": ("rhea2xrefs.tsv", "rhea2xrefs.tsv.gz"),
}

RELEASE_REQUIRED_SOURCES = tuple(SOURCE_BASENAMES)

REACTION_SOURCE_NAMES = (
    "rdf",
    "directions",
    "relationships",
    "obsoletes",
    "reaction_smiles",
)

COMPOUND_SOURCE_NAMES = (
    "sdf",
    "chebi_names",
    "chebi_ph7_3_mapping",
)

CROSS_REFERENCE_SOURCE_NAMES = (
    "xrefs",
    "uniprot_sprot",
    "uniprot_trembl",
)
