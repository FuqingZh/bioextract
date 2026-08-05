from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

_ROOT = Path(__file__).resolve().parents[3]
_THREAD_ENV_VARS = (
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


def test_launcher_clamps_native_environment_before_imports() -> None:
    probe = f"""
import json
import os

import numpy as np
import polars as pl
from threadpoolctl import threadpool_info

np.dot(np.ones((32, 32)), np.ones((32, 32)))
print(json.dumps({{
    "budget": os.environ["BIOEXTRACT_TEST_THREADS"],
    "env": {{name: os.environ[name] for name in {_THREAD_ENV_VARS!r}}},
    "polars": pl.thread_pool_size(),
    "pools": [pool["num_threads"] for pool in threadpool_info()
              if isinstance(pool.get("num_threads"), int)],
}}))
"""
    environment = os.environ.copy()
    environment.pop("BIOEXTRACT_TEST_THREADS", None)
    environment["PYTHONPATH"] = str(_ROOT / "src")
    for name in _THREAD_ENV_VARS:
        environment[name] = "64"

    result = subprocess.run(
        [
            sys.executable,
            str(_ROOT / "scripts/run_with_thread_budget.py"),
            "--default",
            "1",
            "--",
            sys.executable,
            "-c",
            probe,
        ],
        cwd=_ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    payload: dict[str, Any] = json.loads(result.stdout)

    assert payload["budget"] == "1"
    assert payload["env"] == dict.fromkeys(_THREAD_ENV_VARS, "1")
    assert payload["polars"] == 1
    assert all(count <= 1 for count in payload["pools"])
