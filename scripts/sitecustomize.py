"""Set the DuckDB default connection for validation Python interpreters."""

from __future__ import annotations

import os

_BOOTSTRAP_ENABLED = "BIOEXTRACT_RUNTIME_BOOTSTRAP"
_BOOTSTRAP_PID = "BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID"


def _bootstrap_duckdb() -> None:
    current_pid = str(os.getpid())
    if os.environ.get(_BOOTSTRAP_ENABLED) != "1":
        return
    if os.environ.get(_BOOTSTRAP_PID) == current_pid:
        return

    raw_budget = os.environ["BIOEXTRACT_TEST_THREADS"]
    budget = int(raw_budget)
    if budget < 1:
        raise RuntimeError(
            f"BIOEXTRACT_TEST_THREADS must be a positive integer, got {raw_budget!r}"
        )

    import _duckdb

    _duckdb.set_default_connection(_duckdb.connect(config={"threads": str(budget)}))
    os.environ[_BOOTSTRAP_PID] = current_pid


_bootstrap_duckdb()
