#!/usr/bin/env python3
"""Run one repository validation child with its native-thread budget."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "src"))

from bioextract._runtime import (  # noqa: E402
    DEFAULT_TEST_THREADS,
    configure_validation_environment,
)

_BOOTSTRAP_ENABLED = "BIOEXTRACT_RUNTIME_BOOTSTRAP"
_BOOTSTRAP_PID = "BIOEXTRACT_DUCKDB_BOOTSTRAPPED_PID"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--default", type=int, default=DEFAULT_TEST_THREADS)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args()
    if arguments.command and arguments.command[0] == "--":
        arguments.command.pop(0)
    if not arguments.command:
        parser.error("a command is required")
    return arguments


def _is_python_target(command: str) -> bool:
    name = Path(command).name.lower()
    return name in {"pytest", "py.test"} or name.startswith("python")


def _enable_python_bootstrap() -> None:
    scripts_path = str(_ROOT / "scripts")
    python_path = os.environ.get("PYTHONPATH")
    entries = python_path.split(os.pathsep) if python_path else []
    if scripts_path not in entries:
        entries.insert(0, scripts_path)
    os.environ["PYTHONPATH"] = os.pathsep.join(entries)
    os.environ[_BOOTSTRAP_ENABLED] = "1"

    # exec retains this PID but discards the current interpreter and its DuckDB
    # connection. The target interpreter must therefore bootstrap again.
    os.environ.pop(_BOOTSTRAP_PID, None)


def main() -> int:
    arguments = _arguments()
    try:
        configure_validation_environment(default=arguments.default)
    except ValueError as error:
        raise SystemExit(str(error)) from error

    if _is_python_target(arguments.command[0]):
        _enable_python_bootstrap()
    os.execvpe(arguments.command[0], arguments.command, os.environ)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
