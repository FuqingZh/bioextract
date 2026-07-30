"""KEGG metabolic flat-file parsing, publication, and domain access."""

from .core import (
    KEGGMetabolicCapabilityError,
    KEGGMetabolicNamespace,
    KEGGMetabolicSelection,
    MetabolicPublication,
    MetabolicSnapshot,
    evaluate_modules,
    from_metabolic_files,
    from_metabolic_release,
    open_publication,
    validate_selection_namespace,
    write_duckdb,
)

__all__ = [
    "KEGGMetabolicCapabilityError",
    "KEGGMetabolicSelection",
    "KEGGMetabolicNamespace",
    "MetabolicPublication",
    "MetabolicSnapshot",
    "evaluate_modules",
    "from_metabolic_files",
    "from_metabolic_release",
    "open_publication",
    "validate_selection_namespace",
    "write_duckdb",
]
