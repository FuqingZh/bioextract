"""Public database handles exposed through lazy package attributes.

Importing :mod:`bioextract` does not import Polars-backed resource modules.
Each name in ``__all__`` is loaded on first attribute access and then cached in
the package globals. Keep this boundary lazy so callers can inspect or package
``bioextract`` without initializing every resource implementation.
"""

from importlib import import_module
from typing import TYPE_CHECKING, Any

__all__ = [
    "ChEBIDatabase",
    "EggNOGDatabase",
    "GODatabase",
    "InterProDatabase",
    "KEGGDatabase",
    "OmniPathDatabase",
    "ReactomeDatabase",
    "RheaDatabase",
    "STRINGDatabase",
    "UniProtDatabase",
    "WikiPathwaysDatabase",
    "inspect_publication",
]

if TYPE_CHECKING:
    from bioextract.chebi import ChEBIDatabase
    from bioextract.eggnog import EggNOGDatabase
    from bioextract.go import GODatabase
    from bioextract.interpro import InterProDatabase
    from bioextract.kegg import KEGGDatabase
    from bioextract.omnipath import OmniPathDatabase
    from bioextract.publication import inspect_publication
    from bioextract.reactome import ReactomeDatabase
    from bioextract.rhea import RheaDatabase
    from bioextract.stringdb import STRINGDatabase
    from bioextract.uniprot import UniProtDatabase
    from bioextract.wikipathways import WikiPathwaysDatabase

_DATABASE_MODULES: dict[str, str] = {
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

_PUBLIC_MODULES = {**_DATABASE_MODULES, "inspect_publication": "bioextract.publication"}


def __getattr__(name: str) -> Any:
    """Load one declared public object on first attribute access.

    Raises:
        AttributeError: If ``name`` is not declared in ``__all__``.
    """
    module_name = _PUBLIC_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    public_object = getattr(import_module(module_name), name)
    globals()[name] = public_object
    return public_object


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
