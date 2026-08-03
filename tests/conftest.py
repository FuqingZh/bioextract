from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import duckdb
import pytest

_TEST_THREADS_ENV = "BIOEXTRACT_TEST_THREADS"
_DEFAULT_TEST_THREADS = 4


def _test_thread_limit() -> int:
    raw_limit = os.environ.get(_TEST_THREADS_ENV, str(_DEFAULT_TEST_THREADS))
    try:
        limit = int(raw_limit)
    except ValueError as exc:
        message = f"{_TEST_THREADS_ENV} must be a positive integer, got {raw_limit!r}"
        raise pytest.UsageError(message) from exc
    if limit < 1:
        message = f"{_TEST_THREADS_ENV} must be a positive integer, got {raw_limit!r}"
        raise pytest.UsageError(message)
    return limit


_TEST_THREAD_LIMIT = _test_thread_limit()

# Polars reads its thread limit at import time. Root conftest is loaded before
# test modules, so direct `pytest` and repository-owned PDM commands share the
# same bounded default. Callers may still choose a smaller explicit limit.
os.environ.setdefault(_TEST_THREADS_ENV, str(_TEST_THREAD_LIMIT))
os.environ.setdefault("POLARS_MAX_THREADS", str(_TEST_THREAD_LIMIT))
os.environ.setdefault("RAYON_NUM_THREADS", str(_TEST_THREAD_LIMIT))


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
