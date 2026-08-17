from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[3]))

from scripts.check_execution_integrity import check_source, check_tree


def test_current_source_tree_passes_execution_integrity() -> None:
    root = Path(__file__).parents[3] / "src"

    assert check_tree(root) == []


def test_register_io_source_is_restricted_to_shared_adapter() -> None:
    violations = check_source(
        "import polars as pl\npl.io.plugins.register_io_source(source)\n",
        path=Path("src/bioextract/resource.py"),
    )

    assert any("register_io_source" in violation for violation in violations)


def test_unbounded_rows_and_public_limits_are_rejected() -> None:
    violations = check_source(
        "def extract(ids, *, max_rows=None):\n    return list(rows)\n",
        path=Path("src/bioextract/resource.py"),
    )

    assert any("list(rows)" in violation for violation in violations)
    assert any("public method extract" in violation for violation in violations)


def test_bare_type_suppression_is_rejected() -> None:
    violations = check_source(
        "value = unknown()  # type: ignore\n",
        path=Path("src/bioextract/resource.py"),
    )

    assert any("suppression" in violation for violation in violations)
