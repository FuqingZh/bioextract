"""Private runtime configuration for repository-owned validation commands."""

from __future__ import annotations

import os
from collections.abc import MutableMapping
from typing import Final

DEFAULT_TEST_THREADS: Final = 4
NATIVE_THREAD_ENV_VARS: Final[tuple[str, ...]] = (
    "POLARS_MAX_THREADS",
    "RAYON_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMEXPR_MAX_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)
DUCKDB_BOOTSTRAP_PID_ENV: Final = "BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID"


def _positive_integer(name: str, raw_value: str) -> int:
    try:
        value = int(raw_value)
    except ValueError as error:
        raise ValueError(
            f"{name} must be a positive integer, got {raw_value!r}"
        ) from error
    if value < 1:
        raise ValueError(f"{name} must be a positive integer, got {raw_value!r}")
    return value


def configure_validation_environment(
    *,
    default: int = DEFAULT_TEST_THREADS,
    environ: MutableMapping[str, str] | None = None,
) -> int:
    """Validate and apply the native-pool ceiling for one validation process.

    Each unset native-pool variable receives the selected budget. A positive
    caller value is preserved up to that ceiling and clamped above it.

    Returns:
        The selected positive validation budget.

    Raises:
        ValueError: If the default, budget, or any native-pool value is not a
            positive integer.

    Notes:
        Call this before importing analytical libraries. Many native runtimes
        read these variables only while their worker pools are initialized.
    """
    if default < 1:
        raise ValueError(f"default thread budget must be positive, got {default}")

    values = os.environ if environ is None else environ
    budget = _positive_integer(
        "BIOEXTRACT_TEST_THREADS",
        values.get("BIOEXTRACT_TEST_THREADS", str(default)),
    )
    values["BIOEXTRACT_TEST_THREADS"] = str(budget)

    for name in NATIVE_THREAD_ENV_VARS:
        raw_value = values.get(name)
        if raw_value is None:
            values[name] = str(budget)
        else:
            values[name] = str(min(_positive_integer(name, raw_value), budget))

    return budget


def bootstrap_duckdb_default_connection(
    *, environ: MutableMapping[str, str] | None = None
) -> None:
    """Bound DuckDB's default connection once in the current interpreter.

    The process-local marker prevents duplicate setup after ordinary imports
    while allowing spawned children, which have a different PID, to initialize
    their own default connection.
    """
    values = os.environ if environ is None else environ
    current_pid = str(os.getpid())
    if values.get(DUCKDB_BOOTSTRAP_PID_ENV) == current_pid:
        return

    budget = _positive_integer(
        "BIOEXTRACT_TEST_THREADS", values["BIOEXTRACT_TEST_THREADS"]
    )

    import _duckdb

    _duckdb.set_default_connection(_duckdb.connect(config={"threads": str(budget)}))
    values[DUCKDB_BOOTSTRAP_PID_ENV] = current_pid
