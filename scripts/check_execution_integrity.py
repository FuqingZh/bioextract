"""Static checks for the cross-resource execution boundary.

The checker intentionally contains only high-confidence AST and text rules.
It is a repository guard against accidental fake laziness and public execution
knobs; relation-specific streaming behavior remains covered by runtime probes.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

_BARE_IGNORE = re.compile(r"#\s*(?:type:\s*|pyright:\s*)?ignore(?!\[)")
_PUBLIC_LIMIT_NAMES = {"batch_size", "max_rows"}


def check_source(source: str, *, path: Path) -> list[str]:
    """Return integrity violations found in one Python source file."""
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as error:
        return [f"{path}:{error.lineno}: syntax error: {error.msg}"]

    violations: list[str] = []
    adapter_path = path.name == "_lazy.py"
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if (
                isinstance(node.func, ast.Attribute)
                and node.func.attr == "register_io_source"
                and not adapter_path
            ):
                violations.append(
                    f"{path}:{node.lineno}: register_io_source must stay in _lazy.py"
                )
            if (
                isinstance(node.func, ast.Name)
                and node.func.id == "list"
                and len(node.args) == 1
                and isinstance(node.args[0], ast.Name)
                and node.args[0].id == "rows"
            ):
                violations.append(
                    f"{path}:{node.lineno}: list(rows) is an unbounded relation materialization"
                )
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if adapter_path or node.name.startswith("_"):
                continue
            arguments = (
                *node.args.posonlyargs,
                *node.args.args,
                *node.args.kwonlyargs,
            )
            forbidden = sorted(
                argument.arg
                for argument in arguments
                if argument.arg in _PUBLIC_LIMIT_NAMES
            )
            if forbidden:
                violations.append(
                    f"{path}:{node.lineno}: public method {node.name} exposes execution limit(s): "
                    + ", ".join(forbidden)
                )

    for line_number, line in enumerate(source.splitlines(), start=1):
        if "register_deferred_frame_source" in line:
            violations.append(
                f"{path}:{line_number}: removed deferred frame helper reappeared"
            )
        if _BARE_IGNORE.search(line):
            violations.append(
                f"{path}:{line_number}: type suppression must name one diagnostic"
            )
    return violations


def check_tree(root: Path) -> list[str]:
    """Check all Python sources below ``root``."""
    violations: list[str] = []
    for path in sorted(root.rglob("*.py")):
        try:
            source = path.read_text(encoding="utf-8")
        except OSError as error:
            violations.append(f"{path}: cannot read source: {error}")
            continue
        violations.extend(check_source(source, path=path))
    return violations


def main() -> int:
    root = Path(__file__).resolve().parents[1] / "src"
    violations = check_tree(root)
    if violations:
        for violation in violations:
            print(violation, file=sys.stderr)
        return 1
    print(f"execution integrity: checked {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
