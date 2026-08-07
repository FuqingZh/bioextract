"""Set the DuckDB default connection for validation Python interpreters."""

from __future__ import annotations

import os

from bioextract._runtime import bootstrap_duckdb_default_connection

_BOOTSTRAP_ENABLED = "BIOEXTRACT_RUNTIME_BOOTSTRAP"


def _bootstrap_duckdb() -> None:
    if os.environ.get(_BOOTSTRAP_ENABLED) != "1":
        return
    bootstrap_duckdb_default_connection()


_bootstrap_duckdb()
