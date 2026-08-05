"""Private runtime configuration used by repository validation commands."""

from __future__ import annotations

import os
from collections.abc import MutableMapping
from typing import Final

DEFAULT_TEST_THREADS: Final = 4
THREAD_ENV_VARS: Final[tuple[str, ...]] = (
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


def _parse_positive(name: str, raw: str) -> int:
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError(f"{name} must be a positive integer, got {raw!r}") from error
    if value < 1:
        raise ValueError(f"{name} must be a positive integer, got {raw!r}")
    return value


def configure_thread_budget(
    *,
    default: int = DEFAULT_TEST_THREADS,
    environ: MutableMapping[str, str] | None = None,
) -> int:
    """Apply one validation thread ceiling to native-library environments.

    An explicitly smaller native-library value is preserved. Values above the
    common ceiling are clamped so host defaults cannot bypass the repository
    budget. The returned value is the effective ``BIOEXTRACT_TEST_THREADS``
    ceiling.
    """
    if default < 1:
        raise ValueError(f"default thread budget must be positive, got {default}")

    values = os.environ if environ is None else environ
    limit = _parse_positive(
        "BIOEXTRACT_TEST_THREADS",
        values.get("BIOEXTRACT_TEST_THREADS", str(default)),
    )
    values["BIOEXTRACT_TEST_THREADS"] = str(limit)

    for name in THREAD_ENV_VARS:
        raw = values.get(name)
        if raw is None:
            values[name] = str(limit)
            continue
        values[name] = str(min(_parse_positive(name, raw), limit))

    return limit
