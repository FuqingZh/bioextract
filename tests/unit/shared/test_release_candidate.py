from __future__ import annotations

import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[3]))

from scripts.verify_release_candidate import (
    DistributionArtifacts,
    classify_release_tag,
    verify_checkout,
    verify_distribution_versions,
)


def _write_distribution_pair(
    dist: Path,
    tmp_path: Path,
    *,
    wheel_name: str = "bioextract",
    wheel_version: str = "0.8.0",
    sdist_name: str = "bioextract",
    sdist_version: str = "0.8.0",
) -> None:
    dist.mkdir(exist_ok=True)
    wheel = dist / f"bioextract-{wheel_version}-py3-none-any.whl"
    with zipfile.ZipFile(wheel, mode="w") as archive:
        archive.writestr(
            f"bioextract-{wheel_version}.dist-info/METADATA",
            f"Metadata-Version: 2.4\nName: {wheel_name}\nVersion: {wheel_version}\n",
        )
    sdist = dist / f"bioextract-{sdist_version}.tar.gz"
    metadata = tmp_path / f"PKG-INFO-{sdist_version}"
    metadata.write_text(
        f"Metadata-Version: 2.4\nName: {sdist_name}\nVersion: {sdist_version}\n",
        encoding="utf-8",
    )
    with tarfile.open(sdist, mode="w:gz") as archive:
        archive.add(metadata, arcname=f"bioextract-{sdist_version}/PKG-INFO")


@pytest.mark.parametrize("tag", ["0.8.0", "10.2.3"])
def test_final_release_tags_are_canonical(tag: str) -> None:
    assert classify_release_tag(tag, expected_kind="final") == "final"


@pytest.mark.parametrize("tag", ["0.8.0rc1", "10.2.3rc12"])
def test_rc_release_tags_are_canonical(tag: str) -> None:
    assert classify_release_tag(tag, expected_kind="rc") == "rc"


@pytest.mark.parametrize(
    "tag",
    [
        "v0.8.0",
        "0.8",
        "0.8.0rc0",
        "0.8.0a1",
        "0.8.0b1",
        "0.8.0.dev1",
        "0.8.0.post1",
        "0.8.0+local",
        "00.8.0",
    ],
)
def test_other_release_tag_shapes_fail_closed(tag: str) -> None:
    with pytest.raises(ValueError, match="canonical X.Y.Z"):
        classify_release_tag(tag)


def test_github_release_kind_cannot_reclassify_tag() -> None:
    with pytest.raises(ValueError, match="expected final"):
        classify_release_tag("0.8.0rc1", expected_kind="final")
    with pytest.raises(ValueError, match="expected rc"):
        classify_release_tag("0.8.0", expected_kind="rc")


def test_distribution_identity_must_equal_release_tag(tmp_path: Path) -> None:
    dist = tmp_path / "dist"
    _write_distribution_pair(dist, tmp_path)

    artifacts = verify_distribution_versions(dist, expected_version="0.8.0")
    assert artifacts == DistributionArtifacts(
        wheel=dist / "bioextract-0.8.0-py3-none-any.whl",
        sdist=dist / "bioextract-0.8.0.tar.gz",
        version="0.8.0",
    )
    with pytest.raises(ValueError, match="do not match"):
        verify_distribution_versions(dist, expected_version="0.8.1")


def test_distribution_directory_rejects_unverified_artifacts(tmp_path: Path) -> None:
    dist = tmp_path / "dist"
    _write_distribution_pair(dist, tmp_path)
    (dist / "checksums.txt").write_text("unverified\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unexpected artifacts"):
        verify_distribution_versions(dist, expected_version=None)


def test_distribution_metadata_requires_exact_project_name(tmp_path: Path) -> None:
    dist = tmp_path / "dist"
    _write_distribution_pair(dist, tmp_path, wheel_name="not-bioextract")

    with pytest.raises(ValueError, match="Name must be 'bioextract'"):
        verify_distribution_versions(dist, expected_version=None)


def test_wheel_and_sdist_versions_must_agree_without_release_tag(
    tmp_path: Path,
) -> None:
    dist = tmp_path / "dist"
    _write_distribution_pair(dist, tmp_path, sdist_version="0.8.1")

    with pytest.raises(ValueError, match="versions disagree"):
        verify_distribution_versions(dist, expected_version=None)


def test_rc_checkout_must_match_event_tag_and_commit(tmp_path: Path) -> None:
    subprocess.run(["git", "init", "--quiet"], cwd=tmp_path, check=True)
    subprocess.run(
        ["git", "config", "user.email", "tests@bioextract.invalid"],
        cwd=tmp_path,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "bioextract tests"],
        cwd=tmp_path,
        check=True,
    )
    source = tmp_path / "source.txt"
    source.write_text("release candidate\n", encoding="utf-8")
    subprocess.run(["git", "add", "source.txt"], cwd=tmp_path, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "release candidate"],
        cwd=tmp_path,
        check=True,
    )
    subprocess.run(["git", "tag", "0.8.0rc1"], cwd=tmp_path, check=True)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    assert (
        verify_checkout(
            tmp_path,
            tag="0.8.0rc1",
            expected_commit=commit,
        )
        == commit
    )
    with pytest.raises(ValueError, match="event commit"):
        verify_checkout(
            tmp_path,
            tag="0.8.0rc1",
            expected_commit="0" * 40,
        )

    source.write_text("dirty\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must be clean"):
        verify_checkout(
            tmp_path,
            tag="0.8.0rc1",
            expected_commit=commit,
        )
