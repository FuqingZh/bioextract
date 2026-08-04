from __future__ import annotations

import os
from pathlib import Path

import duckdb
import pytest

from bioextract import inspect_publication

_PUBLICATIONS_ENV = "BIOEXTRACT_SMOKE_DUCKDBS"
_PROVENANCE_TABLES = {
    "column_mapping",
    "metadata",
    "source_file",
    "table_info",
    "validation_issue",
}


def _publication_paths() -> tuple[Path, ...]:
    raw_paths = os.environ.get(_PUBLICATIONS_ENV)
    if raw_paths is None:
        pytest.skip(f"{_PUBLICATIONS_ENV} is not configured")
    paths = tuple(
        Path(value).expanduser()
        for value in raw_paths.split(os.pathsep)
        if value.strip()
    )
    if not paths:
        pytest.skip(f"{_PUBLICATIONS_ENV} contains no paths")
    return paths


@pytest.mark.external_snapshot
def test_configured_duckdb_publications_have_embedded_provenance() -> None:
    for path in _publication_paths():
        assert path.is_file(), f"publication is not a file: {path}"
        inspection = inspect_publication(path)

        assert inspection.metadata_schema_version == "1", path
        assert inspection.resource_name, path
        assert inspection.resource_schema_version, path
        assert inspection.source_schema_profile, path
        assert inspection.metadata, path
        assert inspection.tables, path
        assert inspection.table_counts_verified is False, path

        with duckdb.connect(str(path), read_only=True) as connection:
            provenance_tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT table_name FROM information_schema.tables "
                    "WHERE table_schema = '_bioextract'"
                ).fetchall()
            }
        assert provenance_tables == _PROVENANCE_TABLES, path

        print(
            f"{inspection.resource_name}: "
            f"schema={inspection.resource_schema_version} "
            f"release={inspection.release_version} "
            f"status={inspection.validation_status} "
            f"tables={len(inspection.tables)}"
        )
