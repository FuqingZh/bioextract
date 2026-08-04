from __future__ import annotations

import importlib
import subprocess
import sys

import bioextract
import bioextract.errors as errors
from bioextract.errors import BioextractError, CapabilityError, IntegrityError

DATABASE_MODULES = {
    "ChEBIDatabase": "bioextract.chebi",
    "EggNOGDatabase": "bioextract.eggnog",
    "GODatabase": "bioextract.go",
    "InterProDatabase": "bioextract.interpro",
    "KEGGDatabase": "bioextract.kegg",
    "OmniPathDatabase": "bioextract.omnipath",
    "ReactomeDatabase": "bioextract.reactome",
    "RheaDatabase": "bioextract.rhea",
    "STRINGDatabase": "bioextract.stringdb",
    "UniProtDatabase": "bioextract.uniprot",
    "WikiPathwaysDatabase": "bioextract.wikipathways",
}
IMPLEMENTATION_MODULES = {
    class_name: f"{module_name}.{module_name.rsplit('.', maxsplit=1)[-1]}"
    for class_name, module_name in DATABASE_MODULES.items()
}
REMOVED_RESOURCE_EXPORTS = {
    "bioextract.chebi": (
        "ChEBICompoundSelection",
        "ChEBICapabilityError",
        "ChEBIIntegrityError",
    ),
    "bioextract.eggnog": ("EggnogSelection", "EggnogTidyDataset"),
    "bioextract.go": ("GoSubsetId", "GoTidyDataset"),
    "bioextract.interpro": (
        "InterProSelection",
        "InterProTidyConfig",
        "InterProTidyDataset",
    ),
    "bioextract.kegg": ("KeggSelection", "KeggTidyDataset", "KEGGNamespace"),
    "bioextract.omnipath": ("OmniPathSelection",),
    "bioextract.reactome": ("ReactomeSelection", "ReactomeTidyDataset"),
    "bioextract.rhea": (
        "RheaReactionSelection",
        "RheaWriteResult",
        "RheaNamespace",
        "RheaCapabilityError",
    ),
    "bioextract.stringdb": ("StringSelection",),
    "bioextract.uniprot": ("UniProtSelection",),
    "bioextract.wikipathways": (
        "WikiPathwaysSelection",
        "WikiPathwaysTidyDataset",
    ),
}


def test_importing_root_package_does_not_import_resource_modules() -> None:
    module_names = tuple(DATABASE_MODULES.values())
    script = (
        "import sys\n"
        "import bioextract\n"
        f"assert not any(name in sys.modules for name in {module_names!r})\n"
        "assert 'polars' not in sys.modules\n"
        "assert 'duckdb' not in sys.modules\n"
        "assert bioextract.ChEBIDatabase.__name__ == 'ChEBIDatabase'\n"
        "assert 'bioextract.chebi' in sys.modules\n"
        "assert not any(\n"
        f"    name in sys.modules for name in {module_names[1:]!r}\n"
        ")\n"
    )
    subprocess.run([sys.executable, "-c", script], check=True)


def test_top_level_database_exports_are_exact_and_identical() -> None:
    assert bioextract.__all__ == [*DATABASE_MODULES, "inspect_publication"]
    assert set(DATABASE_MODULES) <= set(dir(bioextract))
    for class_name, module_name in DATABASE_MODULES.items():
        resource_module = importlib.import_module(module_name)
        database_class = getattr(resource_module, class_name)
        assert getattr(bioextract, class_name) is database_class
        assert bioextract.__dict__[class_name] is database_class
        assert resource_module.__all__ == [class_name]
        implementation_module = importlib.import_module(
            IMPLEMENTATION_MODULES[class_name]
        )
        assert implementation_module.__all__ == [class_name]


def test_public_errors_are_exact_and_share_one_base_class() -> None:
    assert errors.__all__ == [
        "BioextractError",
        "CapabilityError",
        "IntegrityError",
    ]
    assert issubclass(CapabilityError, BioextractError)
    assert issubclass(IntegrityError, BioextractError)
    assert issubclass(BioextractError, RuntimeError)
    assert not hasattr(bioextract, "BioextractError")
    assert not hasattr(bioextract, "CapabilityError")
    assert not hasattr(bioextract, "IntegrityError")


def test_resource_packages_do_not_reexport_implementation_types() -> None:
    for module_name, removed_names in REMOVED_RESOURCE_EXPORTS.items():
        resource_module = importlib.import_module(module_name)
        for removed_name in removed_names:
            assert not hasattr(resource_module, removed_name)


def test_kegg_implementation_packages_export_no_public_symbols() -> None:
    brite = importlib.import_module("bioextract.kegg.brite")
    metabolic = importlib.import_module("bioextract.kegg.metabolic")
    assert brite.__all__ == []
    assert metabolic.__all__ == []
