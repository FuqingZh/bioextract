from __future__ import annotations

import inspect
import subprocess
import sys

import pytest

import bioextract
import bioextract.publication as publication
from bioextract.errors import IntegrityError


def test_publication_module_has_exact_stable_exports() -> None:
    assert publication.__all__ == [
        "PublicationColumnMapping",
        "PublicationInspection",
        "PublicationMetadata",
        "PublicationSourceFile",
        "PublicationTable",
        "PublicationValidationIssue",
        "inspect_publication",
    ]
    for record_name in publication.__all__[:-1]:
        assert not hasattr(bioextract, record_name)


def test_inspection_is_the_only_additive_lazy_root_export() -> None:
    script = (
        "import sys\n"
        "import bioextract\n"
        "assert 'bioextract.publication' not in sys.modules\n"
        "assert bioextract.inspect_publication.__module__ == "
        "'bioextract.publication'\n"
        "assert 'bioextract.publication' in sys.modules\n"
    )
    subprocess.run([sys.executable, "-c", script], check=True)
    assert bioextract.__all__[-1] == "inspect_publication"


def test_inspection_signature_and_strict_bool_validation() -> None:
    signature = inspect.signature(publication.inspect_publication)
    assert str(signature) == (
        "(path: 'os.PathLike[str] | str', *, verify_table_counts: 'bool' = False) "
        "-> 'PublicationInspection'"
    )
    for invalid in (0, 1, None, "false"):
        with pytest.raises(TypeError, match="must be a bool"):
            publication.inspect_publication(
                "unused.duckdb",
                verify_table_counts=invalid,  # type: ignore[arg-type]
            )


def test_path_resolution_failures_use_public_error_taxonomy() -> None:
    with pytest.raises(IntegrityError, match="invalid.*publication") as caught:
        publication.inspect_publication("invalid\0publication.duckdb")
    assert isinstance(caught.value.__cause__, ValueError)
