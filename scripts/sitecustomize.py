"""Bootstrap bounded DuckDB defaults for validation child interpreters."""

from __future__ import annotations

import os


def _bootstrap_duckdb() -> None:
    if (
        os.environ.get("BIOEXTRACT_RUNTIME_BOOTSTRAP") != "1"
        or os.environ.get("BIOEXTRACT_DUCKDB_BOOTSTRAPPED") == "1"
    ):
        return
    raw_limit = os.environ.get("BIOEXTRACT_TEST_THREADS", "4")
    try:
        limit = int(raw_limit)
    except ValueError as error:
        raise RuntimeError(
            f"BIOEXTRACT_TEST_THREADS must be a positive integer, got {raw_limit!r}"
        ) from error
    if limit < 1:
        raise RuntimeError(
            f"BIOEXTRACT_TEST_THREADS must be a positive integer, got {raw_limit!r}"
        )

    import _duckdb

    _duckdb.set_default_connection(_duckdb.connect(config={"threads": str(limit)}))
    os.environ["BIOEXTRACT_DUCKDB_BOOTSTRAPPED"] = "1"


_bootstrap_duckdb()
