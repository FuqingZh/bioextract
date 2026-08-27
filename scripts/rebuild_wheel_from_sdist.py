from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.verify_release_candidate import (
    DistributionArtifacts,
    verify_distribution_versions,
)


def rebuild_wheel_from_sdist(
    dist: Path,
    *,
    destination: Path,
) -> DistributionArtifacts:
    """Rebuild and verify one wheel from the release sdist without Git metadata."""
    source_artifacts = verify_distribution_versions(dist, expected_version=None)
    with tarfile.open(source_artifacts.sdist, mode="r:gz") as archive:
        if any(".git" in Path(member.name).parts for member in archive.getmembers()):
            raise ValueError("Release sdist must not contain Git metadata")

    if destination.exists():
        if not destination.is_dir() or any(destination.iterdir()):
            raise ValueError(
                "Sdist wheel rebuild destination must be an empty directory"
            )
    else:
        destination.mkdir(parents=True)

    with TemporaryDirectory(prefix="bioextract-uv-state-") as state_root:
        environment = os.environ.copy()
        environment["XDG_CACHE_HOME"] = str(Path(state_root) / "cache")
        environment["UV_CACHE_DIR"] = str(Path(state_root) / "uv-cache")
        completed = subprocess.run(
            [
                "uv",
                "build",
                "--wheel",
                "--no-create-gitignore",
                "--out-dir",
                str(destination),
                str(source_artifacts.sdist),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
            env=environment,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Wheel rebuild from release sdist failed:\n"
            f"{completed.stdout}{completed.stderr}"
        )

    shutil.copy2(source_artifacts.sdist, destination)
    rebuilt_artifacts = verify_distribution_versions(
        destination,
        expected_version=source_artifacts.version,
    )
    print(
        f"rebuilt wheel from sdist without Git metadata: {rebuilt_artifacts.wheel.name}"
    )
    return rebuilt_artifacts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    rebuild_wheel_from_sdist(arguments.dist, destination=arguments.output)


if __name__ == "__main__":
    main()
