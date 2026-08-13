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
from bioextract.chebi._query import ChEBICompoundSelection
from bioextract.eggnog.eggnog import EggnogSelection
from bioextract.go._query import GOAncestorSelection
from bioextract.interpro.interpro import InterProSelection
from bioextract.kegg.kegg import KeggSelection
from bioextract.kegg.metabolic.core import KEGGMetabolicSelection
from bioextract.omnipath.omnipath import OmniPathSelection
from bioextract.reactome.reactome import ReactomeSelection
from bioextract.rhea._query import RheaReactionSelection
from bioextract.stringdb.stringdb import StringSelection
from bioextract.uniprot._query import UniProtSelection
from bioextract.wikipathways.wikipathways import WikiPathwaysSelection


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
    assert not hasattr(kegg.KEGGDatabase, "from_metabolic_release")
    assert not hasattr(chebi.ChEBIDatabase, "from_release")
    assert not hasattr(kegg.KEGGDatabase, "write_parquet")
    assert not hasattr(interpro.InterProDatabase, "write_parquet")


def test_public_relation_classes_do_not_expose_eager_compatibility_methods() -> None:
    relation_classes = (
        chebi.ChEBIDatabase,
        ChEBICompoundSelection,
        eggnog.EggNOGDatabase,
        EggnogSelection,
        go.GODatabase,
        GOAncestorSelection,
        interpro.InterProDatabase,
        InterProSelection,
        kegg.KEGGDatabase,
        KeggSelection,
        KEGGMetabolicSelection,
        omnipath.OmniPathDatabase,
        OmniPathSelection,
        reactome.ReactomeDatabase,
        ReactomeSelection,
        rhea.RheaDatabase,
        RheaReactionSelection,
        stringdb.STRINGDatabase,
        StringSelection,
        uniprot.UniProtDatabase,
        UniProtSelection,
        wikipathways.WikiPathwaysDatabase,
        WikiPathwaysSelection,
    )
    for relation_class in relation_classes:
        public_methods = {
            name
            for name, _member in inspect.getmembers(
                relation_class,
                inspect.isfunction,
            )
            if not name.startswith("_")
        }
        assert not any(name.startswith("extract_") for name in public_methods)
        assert "read_mapping" not in public_methods
        assert "xml_frame" not in public_methods


def test_resource_factories_do_not_expose_limits() -> None:
    factories = (
        go.GODatabase.from_obo,
        kegg.KEGGDatabase.from_brite_json,
        kegg.KEGGDatabase.from_metabolic_files,
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
        chebi.ChEBIDatabase.from_table_files: (
            "source",
            "compounds",
            "names",
            "relations",
            "secondary_ids",
            "database_accessions",
            "structures",
            "chemical_data",
            "chemont_obo",
        ),
        chebi.ChEBIDatabase.from_obo: ("source", "sdf", "chemont_obo"),
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
        kegg.KEGGDatabase.from_metabolic_files: (
            "source",
            "compound_list",
            "compound_entries",
            "reaction_list",
            "reaction_entries",
            "enzyme_list",
            "enzyme_entries",
            "module_list",
            "module_entries",
            "compound_pubchem",
            "compound_reaction",
            "reaction_enzyme",
            "reaction_ko",
            "reaction_module",
            "reaction_pathway",
            "module_pathway",
            "release_version",
        ),
        kegg.KEGGDatabase.from_duckdb: ("path",),
        reactome.ReactomeDatabase.from_files: (
            "uniprot_mapping",
            "pathways",
            "relations",
        ),
        wikipathways.WikiPathwaysDatabase.from_gmt: ("source", "glob"),
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

    metabolic_parameters = inspect.signature(
        kegg.KEGGDatabase.from_metabolic_files
    ).parameters
    assert (
        metabolic_parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    )
    assert metabolic_parameters["source"].default is None
    assert all(
        metabolic_parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        and metabolic_parameters[name].default is None
        for name in tuple(metabolic_parameters)[1:]
    )


def test_phase_2_constructor_parameter_kinds_are_explicit() -> None:
    parameters = inspect.signature(
        wikipathways.WikiPathwaysDatabase.from_gmt
    ).parameters
    assert parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert parameters["glob"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["glob"].default is True


def test_wikipathways_species_scope_is_a_view_operation() -> None:
    parameters = inspect.signature(
        wikipathways.WikiPathwaysDatabase.with_species
    ).parameters
    assert tuple(parameters) == ("self", "species")
    assert parameters["species"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD


def test_go_ancestor_selection_signature_is_explicit() -> None:
    parameters = inspect.signature(go.GODatabase.select_ancestors).parameters
    assert tuple(parameters) == (
        "self",
        "term_ids",
        "target_subset_id",
        "include_self",
        "resolve_alt_ids",
        "include_obsolete",
    )
    assert parameters["term_ids"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert all(
        parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        for name in (
            "target_subset_id",
            "include_self",
            "resolve_alt_ids",
            "include_obsolete",
        )
    )
    assert parameters["target_subset_id"].default is None
    assert parameters["include_self"].default is False
    assert parameters["resolve_alt_ids"].default is True
    assert parameters["include_obsolete"].default is False
