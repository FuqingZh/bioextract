from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PATH_PYPROJECT = ROOT / "pyproject.toml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync the release version into pyproject.toml."
    )
    parser.add_argument(
        "--version",
        required=True,
        help="Release version to write (typically the publish tag).",
    )
    return parser.parse_args()


def replace_once(text: str, pattern: str, replacement: str, context: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Expected to update exactly one {context}, got {count}")
    return updated


def update_pyproject(version: str) -> None:
    text = PATH_PYPROJECT.read_text()
    updated = replace_once(
        text,
        r'(^version = ")[^"]+(")',
        rf"\g<1>{version}\2",
        "project version",
    )
    PATH_PYPROJECT.write_text(updated)


def main() -> None:
    args = parse_args()
    update_pyproject(args.version)
    print(f"Synchronized bioextract release version to {args.version}")


if __name__ == "__main__":
    main()
