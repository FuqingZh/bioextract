from __future__ import annotations

import os
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from bioextract._runtime import configure_validation_environment

try:
    _TEST_THREAD_LIMIT = configure_validation_environment()
except ValueError as error:
    raise pytest.UsageError(str(error)) from error

# Direct pytest does not pass through the repository launcher. Make its
# interpreter and spawned Python children use the same bootstrap boundary.
_SCRIPTS_PATH = str(Path(__file__).resolve().parents[1] / "scripts")
if _SCRIPTS_PATH not in sys.path:
    sys.path.insert(0, _SCRIPTS_PATH)
_python_path = os.environ.get("PYTHONPATH")
_python_path_entries = _python_path.split(os.pathsep) if _python_path else []
if _SCRIPTS_PATH not in _python_path_entries:
    _python_path_entries.insert(0, _SCRIPTS_PATH)
os.environ["PYTHONPATH"] = os.pathsep.join(_python_path_entries)
os.environ["BIOEXTRACT_RUNTIME_BOOTSTRAP"] = "1"

import duckdb  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def bounded_duckdb_threads() -> Iterator[None]:
    original_connect = duckdb.connect

    def connect(
        database: str | Path = ":memory:",
        read_only: bool = False,
        config: dict[str, Any] | None = None,
    ) -> duckdb.DuckDBPyConnection:
        bounded_config = dict(config or {})
        bounded_config.setdefault("threads", str(_TEST_THREAD_LIMIT))
        return original_connect(
            database=database,
            read_only=read_only,
            config=bounded_config,
        )

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(duckdb, "connect", connect)
    yield
    monkeypatch.undo()
