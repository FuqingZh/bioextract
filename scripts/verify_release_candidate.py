from __future__ import annotations

import argparse
import subprocess
import tarfile
import zipfile
from dataclasses import dataclass
from email import message_from_bytes
from email.message import Message
from pathlib import Path
from typing import Literal

from packaging.utils import (
    InvalidSdistFilename,
    InvalidWheelFilename,
    parse_sdist_filename,
    parse_wheel_filename,
)
from packaging.version import InvalidVersion, Version

FINAL_TAG = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
RC_TAG = FINAL_TAG + r"rc[1-9][0-9]*"


@dataclass(frozen=True, slots=True)
class DistributionArtifacts:
    """The only wheel and sdist accepted from one build directory."""

    wheel: Path
    sdist: Path
    version: str


def classify_release_tag(
    tag: str,
    *,
    expected_kind: Literal["auto", "final", "rc"] = "auto",
) -> Literal["final", "rc"]:
    """Validate one canonical project tag and return its release kind."""
    import re

    if re.fullmatch(FINAL_TAG, tag):
        release_kind: Literal["final", "rc"] = "final"
    elif re.fullmatch(RC_TAG, tag):
        release_kind = "rc"
    else:
        raise ValueError(
            "Release tag must be canonical X.Y.Z or X.Y.ZrcN with positive N"
        )

    try:
        parsed = Version(tag)
    except InvalidVersion as error:
        raise ValueError(f"Release tag is not valid PEP 440: {tag!r}") from error
    if str(parsed) != tag:
        raise ValueError(
            f"Release tag must already be canonical: {tag!r} != {str(parsed)!r}"
        )
    if expected_kind != "auto" and release_kind != expected_kind:
        raise ValueError(
            f"Release tag {tag!r} is {release_kind}, expected {expected_kind}"
        )
    return release_kind


def verify_checkout(
    root: Path,
    *,
    tag: str,
    expected_commit: str | None,
) -> str:
    """Assert that a clean checkout is the exact commit referenced by tag."""
    head = _git(root, "rev-parse", "HEAD")
    tag_commit = _git(root, "rev-parse", f"refs/tags/{tag}^{{commit}}")
    if head != tag_commit:
        raise ValueError(
            f"Checkout HEAD {head} does not match release tag {tag!r} at {tag_commit}"
        )
    if expected_commit is not None and expected_commit != tag_commit:
        raise ValueError(
            f"Release event commit {expected_commit} does not match tag commit "
            f"{tag_commit}"
        )
    if _git(root, "status", "--porcelain=v1", "--untracked-files=all"):
        raise ValueError("Release checkout must be clean before building")
    return tag_commit


def verify_distribution_versions(
    dist: Path,
    *,
    expected_version: str | None,
) -> DistributionArtifacts:
    """Validate the complete build directory and both artifact identities."""
    if not dist.is_dir():
        raise ValueError(f"Distribution directory does not exist: {dist}")
    entries = sorted(dist.iterdir())
    if any(not entry.is_file() or entry.is_symlink() for entry in entries):
        raise ValueError("Release dist must contain only regular artifact files")
    wheels = [entry for entry in entries if entry.suffix == ".whl"]
    sdists = [entry for entry in entries if entry.name.endswith(".tar.gz")]
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError("Release build must produce exactly one wheel and one sdist")
    expected_entries = {wheels[0], sdists[0]}
    if set(entries) != expected_entries:
        extras = sorted(entry.name for entry in set(entries) - expected_entries)
        raise ValueError(f"Release dist contains unexpected artifacts: {extras}")

    wheel_name, wheel_version = _wheel_metadata(wheels[0])
    sdist_name, sdist_version = _sdist_metadata(sdists[0])
    observed = {wheels[0].name: wheel_version, sdists[0].name: sdist_version}
    names = {wheels[0].name: wheel_name, sdists[0].name: sdist_name}
    invalid_names = {
        artifact: name for artifact, name in names.items() if name != "bioextract"
    }
    if invalid_names:
        raise ValueError(
            f"Distribution metadata Name must be 'bioextract': {invalid_names}"
        )
    _verify_artifact_filename(wheels[0], name=wheel_name, version=wheel_version)
    _verify_artifact_filename(sdists[0], name=sdist_name, version=sdist_version)
    if len(set(observed.values())) != 1:
        raise ValueError(f"Distribution versions disagree: {observed}")
    mismatched = {
        artifact: version
        for artifact, version in observed.items()
        if expected_version is not None and version != expected_version
    }
    if mismatched:
        raise ValueError(
            f"Distribution versions do not match {expected_version!r}: {mismatched}"
        )
    return DistributionArtifacts(
        wheel=wheels[0],
        sdist=sdists[0],
        version=wheel_version,
    )


def _wheel_metadata(path: Path) -> tuple[str, str]:
    with zipfile.ZipFile(path) as archive:
        metadata_names = [
            name for name in archive.namelist() if name.endswith(".dist-info/METADATA")
        ]
        if len(metadata_names) != 1:
            raise ValueError(f"Wheel must contain exactly one METADATA file: {path}")
        metadata = message_from_bytes(archive.read(metadata_names[0]))
    return _required_identity(metadata, path=path, label="Wheel METADATA")


def _sdist_metadata(path: Path) -> tuple[str, str]:
    with tarfile.open(path, mode="r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.endswith("/PKG-INFO")
        ]
        if len(members) != 1:
            raise ValueError(f"Sdist must contain exactly one PKG-INFO file: {path}")
        handle = archive.extractfile(members[0])
        if handle is None:
            raise ValueError(f"Cannot read sdist PKG-INFO: {path}")
        metadata = message_from_bytes(handle.read())
    return _required_identity(metadata, path=path, label="Sdist PKG-INFO")


def _required_identity(metadata: Message, *, path: Path, label: str) -> tuple[str, str]:
    name = metadata.get("Name")
    version = metadata.get("Version")
    if name is None or version is None:
        raise ValueError(f"{label} has no Name or Version: {path}")
    return str(name), str(version)


def _verify_artifact_filename(path: Path, *, name: str, version: str) -> None:
    try:
        if path.suffix == ".whl":
            filename_name, filename_version, _, _ = parse_wheel_filename(path.name)
        else:
            filename_name, filename_version = parse_sdist_filename(path.name)
    except (InvalidWheelFilename, InvalidSdistFilename) as error:
        raise ValueError(f"Invalid distribution filename: {path.name}") from error
    if filename_name != name or str(filename_version) != version:
        raise ValueError(
            "Distribution filename and embedded metadata disagree: "
            f"{path.name!r} != ({name!r}, {version!r})"
        )


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag")
    parser.add_argument("--kind", choices=("auto", "final", "rc"), default="auto")
    parser.add_argument("--expected-commit")
    parser.add_argument("--dist", type=Path)
    arguments = parser.parse_args()

    if arguments.tag is None:
        if arguments.dist is None:
            parser.error("--tag or --dist is required")
        if arguments.kind != "auto" or arguments.expected_commit is not None:
            parser.error("--kind and --expected-commit require --tag")
        artifacts = verify_distribution_versions(
            arguments.dist,
            expected_version=None,
        )
        print(
            "verified bioextract distribution pair "
            f"{artifacts.version}: {artifacts.wheel.name}, {artifacts.sdist.name}"
        )
        return

    kind = classify_release_tag(arguments.tag, expected_kind=arguments.kind)
    commit = verify_checkout(
        Path.cwd(),
        tag=arguments.tag,
        expected_commit=arguments.expected_commit,
    )
    if arguments.dist is not None:
        verify_distribution_versions(arguments.dist, expected_version=arguments.tag)
    print(f"verified {kind} release candidate {arguments.tag} at {commit}")


if __name__ == "__main__":
    main()
