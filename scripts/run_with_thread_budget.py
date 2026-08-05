"""Run one repository command with the shared native-thread budget."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bioextract._runtime import (
    DEFAULT_TEST_THREADS,
    configure_thread_budget,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--default",
        type=int,
        default=DEFAULT_TEST_THREADS,
        help="budget used when BIOEXTRACT_TEST_THREADS is unset",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command.pop(0)
    if not args.command:
        parser.error("a command is required after the thread-budget options")
    return args


def main() -> int:
    """Configure the environment and replace this process with the command."""
    args = _parse_args()
    try:
        configure_thread_budget(default=args.default)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    command_name = Path(args.command[0]).name.lower()
    if command_name in {"pytest", "py.test"} or command_name.startswith("python"):
        scripts_dir = str(Path(__file__).resolve().parent)
        pythonpath = os.environ.get("PYTHONPATH", "")
        paths = pythonpath.split(os.pathsep) if pythonpath else []
        if scripts_dir not in paths:
            paths.insert(0, scripts_dir)
        os.environ["PYTHONPATH"] = os.pathsep.join(paths)
        os.environ["BIOEXTRACT_RUNTIME_BOOTSTRAP"] = "1"
    os.execvpe(args.command[0], args.command, os.environ)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
