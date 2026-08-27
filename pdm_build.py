from __future__ import annotations

import subprocess
import tomllib
from email import message_from_bytes
from pathlib import Path
from typing import Protocol


class _BuildContext(Protocol):
    root: Path
    target: str


def require_clean_distribution_source(root: Path, *, target: str) -> None:
    """Require clean Git sources or a validated wheel-from-sdist context."""
    if target not in {"wheel", "sdist"}:
        return

    try:
        status = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        if target == "wheel" and _is_valid_sdist_source(root):
            return
        raise RuntimeError(
            "Refusing to build a wheel or sdist because Git cannot identify "
            "one source checkout. Build from a clean Git worktree."
        ) from error
    if status:
        raise RuntimeError(
            "Refusing to build a wheel or sdist from a dirty source tree. "
            "Commit or remove tracked and untracked source changes first."
        )


def _is_valid_sdist_source(root: Path) -> bool:
    """Recognize the static metadata PDM writes into an unpacked sdist."""
    pyproject_path = root / "pyproject.toml"
    package_info_path = root / "PKG-INFO"
    try:
        with pyproject_path.open("rb") as handle:
            project = tomllib.load(handle)["project"]
        metadata = message_from_bytes(package_info_path.read_bytes())
    except (KeyError, OSError, tomllib.TOMLDecodeError):
        return False
    project_name = project.get("name")
    project_version = project.get("version")
    return (
        project_name == "bioextract"
        and isinstance(project_version, str)
        and metadata.get("Name") == project_name
        and metadata.get("Version") == project_version
    )


def pdm_build_initialize(context: _BuildContext) -> None:
    """Reject distribution builds whose source tree is not reproducible."""
    require_clean_distribution_source(context.root, target=context.target)
