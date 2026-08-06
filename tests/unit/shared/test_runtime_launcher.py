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


def test_launcher_bootstraps_duckdb_in_child_python_process() -> None:
    child_probe = """
import json
import os

import duckdb

print(json.dumps({
    "marker": os.environ["BIOEXTRACT_DUCKDB_BOOTSTRAPPED"],
    "pid": os.getpid(),
    "threads": duckdb.sql("SELECT current_setting('threads')").fetchone()[0],
}))
"""
    parent_probe = f"""
import subprocess
import sys

subprocess.run([sys.executable, "-c", {child_probe!r}], check=True)
    """
    environment = os.environ.copy()
    environment.pop("BIOEXTRACT_TEST_THREADS", None)
    environment["PYTHONPATH"] = str(_ROOT / "src")
    result = subprocess.run(
        [
            sys.executable,
            str(_ROOT / "scripts/run_with_thread_budget.py"),
            "--default",
            "1",
            "--",
            sys.executable,
            "-c",
            parent_probe,
        ],
        cwd=_ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    payload: dict[str, Any] = json.loads(result.stdout)

    assert payload == {
        "marker": str(payload["pid"]),
        "pid": payload["pid"],
        "threads": 1,
    }


def test_launcher_rebootstraps_duckdb_after_same_pid_exec() -> None:
    probe = """
import json
import os

import duckdb

print(json.dumps({
    "marker": os.environ["BIOEXTRACT_DUCKDB_BOOTSTRAPPED"],
    "pid": os.getpid(),
    "threads": duckdb.sql("SELECT current_setting('threads')").fetchone()[0],
}))
"""
    environment = os.environ.copy()
    environment.pop("BIOEXTRACT_TEST_THREADS", None)
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(_ROOT / "scripts"), str(_ROOT / "src"))
    )
    environment["BIOEXTRACT_RUNTIME_BOOTSTRAP"] = "1"
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

    assert payload == {
        "marker": str(payload["pid"]),
        "pid": payload["pid"],
        "threads": 1,
    }
