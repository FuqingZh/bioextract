"""Public resource namespaces exposed through lazy package attributes.

Importing :mod:`bioextract` does not import Polars-backed resource modules.
Each name in ``__all__`` is loaded on first attribute access and then cached in
the package globals. Keep this boundary lazy so callers can inspect or package
``bioextract`` without initializing every resource implementation.
"""

from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any

__all__ = [
    "chebi",
    "eggnog",
    "go",
    "interpro",
    "kegg",
    "omnipath",
    "reactome",
    "rhea",
    "stringdb",
    "uniprot",
    "wikipathways",
]

if TYPE_CHECKING:
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

_ALIAS_MODULES: dict[str, str] = {
    "chebi": "bioextract.chebi",
    "eggnog": "bioextract.eggnog",
    "go": "bioextract.go",
    "interpro": "bioextract.interpro",
    "kegg": "bioextract.kegg",
    "omnipath": "bioextract.omnipath",
    "reactome": "bioextract.reactome",
    "rhea": "bioextract.rhea",
    "stringdb": "bioextract.stringdb",
    "uniprot": "bioextract.uniprot",
    "wikipathways": "bioextract.wikipathways",
}


def __getattr__(name: str) -> Any:
    """Load one declared resource namespace on first attribute access.

    Raises:
        AttributeError: If ``name`` is not one of the namespaces in ``__all__``.
    """
    module_name = _ALIAS_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_loaded: ModuleType = import_module(module_name)
    globals()[name] = module_loaded
    return module_loaded


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
