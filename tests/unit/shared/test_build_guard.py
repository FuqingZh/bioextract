from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tarfile
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from packaging.version import Version

sys.path.insert(0, str(Path(__file__).parents[3]))

from scripts.build_guard import require_clean_distribution_source
from scripts.verify_release_candidate import verify_distribution_versions

REPOSITORY_ROOT = Path(__file__).parents[3]


def _initialize_repository(path: Path) -> None:
    subprocess.run(["git", "init", "--quiet"], cwd=path, check=True)
    subprocess.run(
        ["git", "config", "user.email", "tests@bioextract.invalid"],
        cwd=path,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "bioextract tests"],
        cwd=path,
        check=True,
    )
    (path / "tracked.txt").write_text("clean\n", encoding="utf-8")
    subprocess.run(["git", "add", "tracked.txt"], cwd=path, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "test fixture"],
        cwd=path,
        check=True,
    )


@pytest.mark.parametrize("target", ["wheel", "sdist"])
@pytest.mark.parametrize("dirty_kind", ["tracked", "untracked"])
def test_distribution_build_rejects_dirty_source_before_artifact_creation(
    tmp_path: Path,
    target: str,
    dirty_kind: str,
) -> None:
    _initialize_repository(tmp_path)
    if dirty_kind == "tracked":
        (tmp_path / "tracked.txt").write_text("dirty\n", encoding="utf-8")
    else:
        (tmp_path / "untracked.txt").write_text("dirty\n", encoding="utf-8")
    dist = tmp_path / "dist"
    dist.mkdir()

    with pytest.raises(RuntimeError, match="dirty source tree"):
        require_clean_distribution_source(tmp_path, target=target)

    assert not list(dist.iterdir())


def test_clean_distribution_build_and_dirty_editable_are_allowed(
    tmp_path: Path,
) -> None:
    _initialize_repository(tmp_path)
    require_clean_distribution_source(tmp_path, target="wheel")

    (tmp_path / "tracked.txt").write_text("dirty\n", encoding="utf-8")
    require_clean_distribution_source(tmp_path, target="editable")


def test_distribution_guard_rejects_source_without_scm_checkout(
    tmp_path: Path,
) -> None:
    with pytest.raises(RuntimeError, match="Git cannot identify one source checkout"):
        require_clean_distribution_source(tmp_path, target="sdist")


def _write_actual_build_project(path: Path) -> None:
    package = path / "src" / "bioextract"
    scripts = path / "scripts"
    package.mkdir(parents=True)
    scripts.mkdir()
    package.joinpath("__init__.py").write_text(
        '"""Actual PDM build fixture."""\n',
        encoding="utf-8",
    )
    path.joinpath("pyproject.toml").write_text(
        """[project]
name = "bioextract"
dynamic = ["version"]
description = "PDM build integration fixture"
requires-python = ">=3.13"

[build-system]
requires = ["pdm-backend>=2.4"]
build-backend = "pdm.backend"

[tool.pdm]
distribution = true

[tool.pdm.version]
source = "scm"
tag_filter = "*.*.*"
tag_regex = '^(?P<version>(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)(?:rc[1-9][0-9]*)?)$'
""",
        encoding="utf-8",
    )
    path.joinpath(".gitignore").write_text(
        "/dist/\n/.pdm-build/\n*.egg-info/\n__pycache__/\n",
        encoding="utf-8",
    )
    path.joinpath("pdm_build.py").write_text(
        REPOSITORY_ROOT.joinpath("pdm_build.py").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    scripts.joinpath("build_guard.py").write_text(
        REPOSITORY_ROOT.joinpath("scripts/build_guard.py").read_text(encoding="utf-8"),
        encoding="utf-8",
    )


def _commit_build_project(path: Path, *, tag: str | None = None) -> None:
    subprocess.run(["git", "init", "--quiet"], cwd=path, check=True)
    subprocess.run(
        ["git", "config", "user.email", "tests@bioextract.invalid"],
        cwd=path,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "bioextract tests"],
        cwd=path,
        check=True,
    )
    subprocess.run(["git", "add", "."], cwd=path, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "build fixture"],
        cwd=path,
        check=True,
    )
    if tag is not None:
        subprocess.run(["git", "tag", tag], cwd=path, check=True)


def _actual_pdm_build(
    path: Path,
    *,
    extra_arguments: Sequence[str] = (),
) -> subprocess.CompletedProcess[str]:
    with TemporaryDirectory(prefix="bioextract-pdm-state-") as state_root:
        environment = os.environ.copy()
        environment["XDG_STATE_HOME"] = state_root
        environment["XDG_CACHE_HOME"] = str(Path(state_root) / "cache")
        environment["UV_CACHE_DIR"] = str(Path(state_root) / "uv-cache")
        return subprocess.run(
            ["pdm", "build", "--dest", "dist", *extra_arguments],
            cwd=path,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
            env=environment,
        )


def _actual_wheel_build_from_sdist(
    sdist: Path,
    *,
    destination: Path,
) -> subprocess.CompletedProcess[str]:
    with TemporaryDirectory(prefix="bioextract-uv-state-") as state_root:
        environment = os.environ.copy()
        environment["XDG_CACHE_HOME"] = str(Path(state_root) / "cache")
        environment["UV_CACHE_DIR"] = str(Path(state_root) / "uv-cache")
        return subprocess.run(
            [
                "uv",
                "build",
                "--wheel",
                "--no-create-gitignore",
                "--out-dir",
                str(destination),
                str(sdist),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
            env=environment,
        )


@pytest.mark.parametrize("tag", ["0.8.0", "0.8.0rc1"])
def test_actual_isolated_pdm_build_uses_scm_tag_and_project_hook(
    tmp_path: Path,
    tag: str,
) -> None:
    _write_actual_build_project(tmp_path)
    _commit_build_project(tmp_path, tag=tag)

    completed = _actual_pdm_build(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    artifacts = verify_distribution_versions(
        tmp_path / "dist",
        expected_version=tag,
    )
    assert artifacts.version == tag


def test_actual_untagged_pdm_build_has_consistent_scm_identity(tmp_path: Path) -> None:
    _write_actual_build_project(tmp_path)
    _commit_build_project(tmp_path)

    completed = _actual_pdm_build(tmp_path, extra_arguments=("--no-isolation",))

    assert completed.returncode == 0, completed.stdout + completed.stderr
    artifacts = verify_distribution_versions(
        tmp_path / "dist",
        expected_version=None,
    )
    assert Version(artifacts.version).is_devrelease is True


def test_actual_pdm_wheel_rebuilds_from_sdist_without_git(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _write_actual_build_project(source)
    _commit_build_project(source, tag="0.8.0rc1")

    source_build = _actual_pdm_build(source)

    assert source_build.returncode == 0, source_build.stdout + source_build.stderr
    source_artifacts = verify_distribution_versions(
        source / "dist",
        expected_version="0.8.0rc1",
    )

    with tarfile.open(source_artifacts.sdist, mode="r:gz") as archive:
        assert all(".git" not in Path(member.name).parts for member in archive)

    rebuilt_dist = tmp_path / "rebuilt-dist"
    wheel_build = _actual_wheel_build_from_sdist(
        source_artifacts.sdist,
        destination=rebuilt_dist,
    )

    assert wheel_build.returncode == 0, wheel_build.stdout + wheel_build.stderr
    shutil.copy2(source_artifacts.sdist, rebuilt_dist)
    rebuilt_artifacts = verify_distribution_versions(
        rebuilt_dist,
        expected_version="0.8.0rc1",
    )
    assert rebuilt_artifacts.version == source_artifacts.version


@pytest.mark.parametrize("dirty_kind", ["tracked", "untracked"])
def test_actual_pdm_build_rejects_dirty_checkout(
    tmp_path: Path,
    dirty_kind: str,
) -> None:
    _write_actual_build_project(tmp_path)
    _commit_build_project(tmp_path, tag="0.8.0")
    if dirty_kind == "tracked":
        (tmp_path / "src" / "bioextract" / "__init__.py").write_text(
            '"""Dirty fixture."""\n',
            encoding="utf-8",
        )
    else:
        (tmp_path / "untracked.py").write_text("dirty = True\n", encoding="utf-8")

    completed = _actual_pdm_build(tmp_path, extra_arguments=("--no-isolation",))

    assert completed.returncode != 0
    assert "dirty source tree" in completed.stdout + completed.stderr
    assert not list((tmp_path / "dist").glob("*.whl"))
    assert not list((tmp_path / "dist").glob("*.tar.gz"))


def test_actual_pdm_build_rejects_source_without_scm_checkout(tmp_path: Path) -> None:
    _write_actual_build_project(tmp_path)

    completed = _actual_pdm_build(tmp_path, extra_arguments=("--no-isolation",))

    assert completed.returncode != 0
    assert not list((tmp_path / "dist").glob("*.whl"))
    assert not list((tmp_path / "dist").glob("*.tar.gz"))
