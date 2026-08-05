from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from bioextract._runtime import configure_thread_budget

try:
    _TEST_THREAD_LIMIT = configure_thread_budget()
except ValueError as error:
    raise pytest.UsageError(str(error)) from error

# Make direct pytest subprocesses inherit the same early bootstrap as the PDM
# launcher. The current interpreter still needs the low-level setup when it
# was not started through that launcher.
os.environ["BIOEXTRACT_RUNTIME_BOOTSTRAP"] = "1"
_SCRIPTS_DIR = str(Path(__file__).resolve().parents[1] / "scripts")
_pythonpath = os.environ.get("PYTHONPATH", "")
_pythonpath_parts = _pythonpath.split(os.pathsep) if _pythonpath else []
if _SCRIPTS_DIR not in _pythonpath_parts:
    _pythonpath_parts.insert(0, _SCRIPTS_DIR)
os.environ["PYTHONPATH"] = os.pathsep.join(_pythonpath_parts)

if os.environ.get("BIOEXTRACT_DUCKDB_BOOTSTRAPPED") != "1":
    # Import the low-level extension first so its bounded connection becomes
    # the package default. The high-level package otherwise creates an
    # unbounded default scheduler before the connection wrapper runs.
    import _duckdb

    _duckdb.set_default_connection(
        _duckdb.connect(config={"threads": str(_TEST_THREAD_LIMIT)})
    )

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
