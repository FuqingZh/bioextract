from __future__ import annotations

import os
from pathlib import Path

import duckdb
import pytest

_PUBLICATIONS_ENV = "BIOEXTRACT_SMOKE_DUCKDBS"


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
        with duckdb.connect(str(path), read_only=True) as connection:
            internal_schema = connection.execute(
                "SELECT count(*) FROM information_schema.schemata "
                "WHERE schema_name = '_bioextract'"
            ).fetchone()
            metadata_rows = connection.execute(
                "SELECT count(*) FROM _bioextract.metadata"
            ).fetchone()
        assert internal_schema == (1,), path
        assert metadata_rows is not None
        assert metadata_rows[0] > 0, path
