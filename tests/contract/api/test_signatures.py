from __future__ import annotations

import inspect

import pytest

import bioextract.chebi as chebi
import bioextract.eggnog as eggnog
import bioextract.go as go
import bioextract.interpro as interpro
import bioextract.kegg as kegg
import bioextract.omnipath as omnipath
import bioextract.reactome as reactome
import bioextract.rhea as rhea
import bioextract.stringdb as stringdb
import bioextract.uniprot as uniprot
import bioextract.wikipathways as wikipathways
from bioextract._tidy import TidyDataset
from bioextract.stringdb.stringdb import StringSelection


@pytest.mark.parametrize(
    ("module", "old_name"),
    [
        (go, "GoDb"),
        (kegg, "KeggDb"),
        (reactome, "ReactomeDb"),
        (wikipathways, "WikiPathwaysDb"),
        (eggnog, "EggnogDb"),
        (interpro, "InterProDb"),
        (uniprot, "UniprotDb"),
        (stringdb, "StringDb"),
        (omnipath, "OmniPathDb"),
        (rhea, "RheaDb"),
    ],
)
def test_legacy_database_type_aliases_are_not_exported(
    module: object,
    old_name: str,
) -> None:
    assert not hasattr(module, old_name)


def test_removed_legacy_writer_apis_are_not_exported() -> None:
    assert not hasattr(TidyDataset, "write")
    assert not hasattr(TidyDataset, "write_parquet")
    assert not hasattr(StringSelection, "with_score_min")
    assert not hasattr(eggnog.EggNOGDatabase, "from_files")
    assert not hasattr(eggnog.EggNOGDatabase, "extract_mapping")
    assert not hasattr(eggnog.EggNOGDatabase, "write_parquet")
    assert not hasattr(eggnog.EggNOGDatabase, "write_duckdb")
    assert not hasattr(eggnog.EggNOGDatabase, "write_tidy")
    assert not hasattr(eggnog.EggNOGDatabase, "from_duckdb")
    assert not hasattr(eggnog.EggNOGDatabase, "persist")
    assert not hasattr(eggnog.EggNOGDatabase, "unpack")
    assert not hasattr(kegg.KEGGDatabase, "write_parquet")
    assert not hasattr(interpro.InterProDatabase, "write_parquet")


def test_resource_factories_do_not_expose_limits() -> None:
    factories = (
        chebi.ChEBIDatabase.from_release,
        go.GODatabase.from_obo,
        kegg.KEGGDatabase.from_brite_json,
        reactome.ReactomeDatabase.from_files,
        wikipathways.WikiPathwaysDatabase.from_gmt,
        eggnog.EggNOGDatabase.from_sqlite,
        interpro.InterProDatabase.from_mapping_files,
        uniprot.UniProtDatabase.from_idmapping,
        uniprot.UniProtDatabase.from_knowledgebase,
        stringdb.STRINGDatabase.from_files,
        omnipath.OmniPathDatabase.from_files,
        rhea.RheaDatabase.from_files,
    )
    assert all(
        "limits" not in inspect.signature(factory).parameters for factory in factories
    )


def test_resource_factory_parameter_names_follow_domain_roles() -> None:
    expected = {
        chebi.ChEBIDatabase.from_release: ("source", "chemont_obo"),
        chebi.ChEBIDatabase.from_table_files: (
            "compounds",
            "names",
            "relations",
            "secondary_ids",
            "database_accessions",
            "structures",
            "chemical_data",
            "chemont_obo",
        ),
        chebi.ChEBIDatabase.from_obo: ("path", "sdf", "chemont_obo"),
        chebi.ChEBIDatabase.from_duckdb: ("path",),
        go.GODatabase.from_obo: ("path",),
        go.GODatabase.from_duckdb: ("path",),
        kegg.KEGGDatabase.from_brite_json: ("path",),
        kegg.KEGGDatabase.from_mapping_files: (
            "uniprot_conversion",
            "gene_ko",
            "gene_pathway",
            "organism_code",
            "gene_list",
            "ncbi_gene_conversion",
        ),
        kegg.KEGGDatabase.from_metabolic_release: ("source", "release_version"),
        kegg.KEGGDatabase.from_duckdb: ("path",),
        reactome.ReactomeDatabase.from_files: (
            "uniprot_mapping",
            "pathways",
            "relations",
        ),
        wikipathways.WikiPathwaysDatabase.from_gmt: ("source", "species", "glob"),
        wikipathways.WikiPathwaysDatabase.from_duckdb: ("path",),
        eggnog.EggNOGDatabase.from_sqlite: (
            "source",
            "cog_functions",
            "temp_dir",
        ),
        interpro.InterProDatabase.from_mapping_files: (
            "protein_to_interpro",
            "interpro_xml",
        ),
        interpro.InterProDatabase.from_duckdb: ("path",),
        uniprot.UniProtDatabase.from_idmapping: ("path", "release_version"),
        uniprot.UniProtDatabase.from_knowledgebase: (
            "entries",
            "canonical_sequences",
            "isoform_sequences",
            "release_version",
        ),
        uniprot.UniProtDatabase.from_duckdb: ("path",),
        stringdb.STRINGDatabase.from_files: (
            "aliases",
            "links",
            "rank_by_source",
            "release_version",
        ),
        omnipath.OmniPathDatabase.from_files: ("enzsub", "interactions"),
        rhea.RheaDatabase.from_files: (
            "source",
            "rdf",
            "directions",
            "relationships",
            "obsolete_reactions",
            "reaction_smiles",
            "sdf",
            "chebi_names",
            "chebi_ph7_3_mapping",
            "xrefs",
            "uniprot_sprot",
            "uniprot_trembl",
        ),
        rhea.RheaDatabase.from_duckdb: ("path",),
    }
    assert {
        factory: tuple(inspect.signature(factory).parameters) for factory in expected
    } == expected


def test_phase_1_constructor_parameter_kinds_are_explicit() -> None:
    eggnog_parameters = inspect.signature(eggnog.EggNOGDatabase.from_sqlite).parameters
    assert eggnog_parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert all(
        eggnog_parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        for name in ("cog_functions", "temp_dir")
    )
    assert "cache_dir" not in eggnog_parameters
    assert "persist" not in eggnog_parameters

    string_parameters = inspect.signature(stringdb.STRINGDatabase.from_files).parameters
    assert all(
        parameter.kind is inspect.Parameter.KEYWORD_ONLY
        for parameter in string_parameters.values()
    )

    rhea_parameters = inspect.signature(rhea.RheaDatabase.from_files).parameters
    assert rhea_parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert rhea_parameters["source"].default is None
    assert all(
        rhea_parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        and rhea_parameters[name].default is None
        for name in tuple(rhea_parameters)[1:]
    )


def test_phase_2_constructor_parameter_kinds_are_explicit() -> None:
    parameters = inspect.signature(
        wikipathways.WikiPathwaysDatabase.from_gmt
    ).parameters
    assert parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert all(
        parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        for name in ("species", "glob")
    )
    assert parameters["glob"].default is True
