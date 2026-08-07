from __future__ import annotations

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any, cast

import pytest

from bioextract._runtime import NATIVE_THREAD_ENV_VARS

_ROOT = Path(__file__).resolve().parents[3]
_LAUNCHER = _ROOT / "scripts" / "run_with_thread_budget.py"
_SCRIPTS = _ROOT / "scripts"


def _clean_environment() -> dict[str, str]:
    environment = os.environ.copy()
    for name in (
        "BIOEXTRACT_TEST_THREADS",
        "BIOEXTRACT_RUNTIME_BOOTSTRAP",
        "BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID",
        *NATIVE_THREAD_ENV_VARS,
    ):
        environment.pop(name, None)
    environment["PYTHONPATH"] = str(_ROOT / "src")
    return environment


def _launch_python(
    source: str,
    *,
    environment: dict[str, str] | None = None,
    default: int = 4,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(_LAUNCHER),
            "--default",
            str(default),
            "--",
            sys.executable,
            "-c",
            source,
        ],
        cwd=_ROOT,
        env=environment or _clean_environment(),
        check=True,
        capture_output=True,
        text=True,
    )


def test_launcher_applies_environment_to_non_python_child() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(_LAUNCHER),
            "--default",
            "2",
            "--",
            "/usr/bin/env",
        ],
        cwd=_ROOT,
        env=_clean_environment(),
        check=True,
        capture_output=True,
        text=True,
    )
    child_environment = dict(
        line.split("=", 1) for line in result.stdout.splitlines() if "=" in line
    )

    assert child_environment["BIOEXTRACT_TEST_THREADS"] == "2"
    assert all(child_environment[name] == "2" for name in NATIVE_THREAD_ENV_VARS)
    assert "BIOEXTRACT_RUNTIME_BOOTSTRAP" not in child_environment


def test_repository_validation_children_use_budget_launcher() -> None:
    configuration = tomllib.loads((_ROOT / "pyproject.toml").read_text())
    scripts = cast(
        dict[str, str | dict[str, str]],
        configuration["tool"]["pdm"]["scripts"],
    )

    for name in (
        "lock-check",
        "format-check",
        "lint",
        "typecheck",
        "test-unit",
        "test-contract",
        "test-integration",
        "test-smoke",
        "test",
    ):
        definition = scripts[name]
        command = definition["cmd"] if isinstance(definition, dict) else definition
        assert command.startswith("python scripts/run_with_thread_budget.py "), name


def test_launcher_bootstraps_spawned_python_child() -> None:
    child_probe = """
import json
import os
import duckdb
print(json.dumps({
    "marker": os.environ["BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID"],
    "pid": os.getpid(),
    "threads": duckdb.sql("SELECT current_setting('threads')").fetchone()[0],
}))
"""
    parent_probe = f"""
import subprocess
import sys
subprocess.run([sys.executable, "-c", {child_probe!r}], check=True)
"""

    payload = json.loads(_launch_python(parent_probe, default=1).stdout)

    assert payload["marker"] == str(payload["pid"])
    assert payload["threads"] == 1


def test_launcher_clears_same_pid_marker_before_python_exec() -> None:
    probe = """
import json
import os
import duckdb
print(json.dumps({
    "marker": os.environ["BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID"],
    "pid": os.getpid(),
    "threads": duckdb.sql("SELECT current_setting('threads')").fetchone()[0],
}))
"""
    environment = _clean_environment()
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(_SCRIPTS), environment["PYTHONPATH"])
    )
    environment["BIOEXTRACT_RUNTIME_BOOTSTRAP"] = "1"

    payload = json.loads(
        _launch_python(probe, environment=environment, default=1).stdout
    )

    assert payload["marker"] == str(payload["pid"])
    assert payload["threads"] == 1


def test_launcher_bounds_loaded_native_pools() -> None:
    probe = f"""
import json
import os
import duckdb
import numpy as np
import polars as pl
from threadpoolctl import threadpool_info

np.dot(np.ones((32, 32)), np.ones((32, 32)))
duckdb_threads = duckdb.sql(
    "SELECT current_setting('threads')::INTEGER"
).fetchone()[0]
print(json.dumps({{
    "environment": {{name: os.environ[name] for name in {NATIVE_THREAD_ENV_VARS!r}}},
    "polars": pl.thread_pool_size(),
    "duckdb": duckdb_threads,
    "pools": threadpool_info(),
}}))
"""

    payload: dict[str, Any] = json.loads(_launch_python(probe, default=1).stdout)

    assert payload["environment"] == dict.fromkeys(NATIVE_THREAD_ENV_VARS, "1")
    assert payload["polars"] <= 1
    assert payload["duckdb"] <= 1
    for pool in payload["pools"]:
        assert pool["num_threads"] <= 1, pool


@pytest.mark.parametrize("invalid", ["0", "-1", "many"])
def test_launcher_rejects_invalid_budget_before_child(invalid: str) -> None:
    environment = _clean_environment()
    environment["BIOEXTRACT_TEST_THREADS"] = invalid

    result = subprocess.run(
        [
            sys.executable,
            str(_LAUNCHER),
            "--",
            sys.executable,
            "-c",
            "print('child started')",
        ],
        cwd=_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "must be a positive integer" in result.stderr
    assert "child started" not in result.stdout
