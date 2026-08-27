from __future__ import annotations

from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory

from bioextract import ChEBIDatabase, inspect_publication

_MINIMAL_CHEBI_OBO = """format-version: 1.2

[Term]
id: CHEBI:1
name: water
"""


def main() -> None:
    installed_version = version("bioextract")
    with TemporaryDirectory(prefix="bioextract-wheel-smoke-") as temp_directory:
        root = Path(temp_directory)
        source = root / "chebi.obo"
        source.write_text(_MINIMAL_CHEBI_OBO, encoding="utf-8")
        publication_path = root / "chebi.duckdb"
        ChEBIDatabase.from_obo(source).write_duckdb(publication_path)
        inspected = inspect_publication(publication_path)
    if inspected.package_version != installed_version:
        raise RuntimeError(
            "Installed wheel metadata and publication metadata differ: "
            f"{installed_version!r} != {inspected.package_version!r}"
        )
    print(f"installed wheel publication identity verified: {installed_version}")


if __name__ == "__main__":
    main()
