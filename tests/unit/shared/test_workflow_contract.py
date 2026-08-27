from __future__ import annotations

import re
from pathlib import Path
from typing import cast

import pytest
import yaml

REPOSITORY_ROOT = Path(__file__).parents[3]
WORKFLOW_ROOT = REPOSITORY_ROOT / ".github" / "workflows"
PYPI_ACTION = "pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33"
UPLOAD_ARTIFACT_ACTION = (
    "actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a"
)
DOWNLOAD_ARTIFACT_ACTION = (
    "actions/download-artifact@3e5f45b2cfb9172054b4087a40e8e0b5a5461e7c"
)


def _workflow(name: str) -> str:
    return (WORKFLOW_ROOT / name).read_text(encoding="utf-8")


def _mapping(value: object) -> dict[str, object]:
    assert isinstance(value, dict)
    return cast("dict[str, object]", value)


def _workflow_document(name: str) -> dict[str, object]:
    document = yaml.load(_workflow(name), Loader=yaml.BaseLoader)
    return _mapping(document)


def _steps(job: object) -> list[dict[str, object]]:
    raw_steps = _mapping(job)["steps"]
    assert isinstance(raw_steps, list)
    return [_mapping(step) for step in cast("list[object]", raw_steps)]


@pytest.mark.parametrize(
    "name",
    ["py-ci.yml", "nightly.yml", "release-build.yml", "publish.yml"],
)
def test_workflows_are_structurally_valid_yaml(name: str) -> None:
    document = _workflow_document(name)
    assert isinstance(document["name"], str)
    assert _mapping(document["on"])
    assert _mapping(document["permissions"]) == {"contents": "read"}
    assert _mapping(document["jobs"])


@pytest.mark.parametrize(
    "name",
    ["py-ci.yml", "nightly.yml", "release-build.yml", "publish.yml"],
)
def test_workflow_actions_are_pinned_to_full_commit_sha(name: str) -> None:
    jobs = _mapping(_workflow_document(name)["jobs"])
    for job in jobs.values():
        for step in _steps(job):
            uses = step.get("uses")
            if uses is None:
                continue
            assert re.fullmatch(r"[^@]+@[0-9a-f]{40}", cast("str", uses))


def test_only_published_github_releases_can_upload_to_pypi() -> None:
    publish = _workflow("publish.yml")
    document = _workflow_document("publish.yml")
    events = _mapping(document["on"])
    assert set(events) == {"release"}
    assert _mapping(events["release"])["types"] == ["published"]
    jobs = _mapping(document["jobs"])
    assert set(jobs) == {
        "final_build",
        "final_publish",
        "rc_build",
        "rc_publish",
    }
    for release_kind, prerelease_expression in (
        ("final", "!github.event.release.prerelease"),
        ("rc", "github.event.release.prerelease"),
    ):
        build_job_name = f"{release_kind}_build"
        build_job = _mapping(jobs[build_job_name])
        assert prerelease_expression in cast("str", build_job["if"])
        assert _mapping(build_job["permissions"]) == {"contents": "read"}
        assert "environment" not in build_job
        build_steps = _steps(build_job)
        assert all(step.get("uses") != PYPI_ACTION for step in build_steps)
        artifact_steps = [
            step for step in build_steps if step.get("uses") == UPLOAD_ARTIFACT_ACTION
        ]
        assert len(artifact_steps) == 1
        artifact_with = _mapping(artifact_steps[0]["with"])
        assert artifact_with["path"] == "dist/*"
        assert "github.event.release.tag_name" in cast("str", artifact_with["name"])
        verification_steps = [
            step
            for step in build_steps
            if str(step.get("name", "")).startswith("Verify ")
        ]
        assert any(
            "--dist dist" in str(step.get("run", "")) for step in verification_steps
        )
        assert any(
            step.get("name") == "Rebuild wheel from sdist"
            and "scripts.rebuild_wheel_from_sdist" in str(step.get("run", ""))
            for step in build_steps
        )
        assert any(
            "build/sdist-wheel/*.whl" in str(step.get("run", ""))
            for step in build_steps
        )

        publish_job = _mapping(jobs[f"{release_kind}_publish"])
        publish_condition = cast("str", publish_job["if"])
        assert prerelease_expression in publish_condition
        assert f"needs.{build_job_name}.result == 'success'" in publish_condition
        assert publish_job["needs"] == build_job_name
        assert _mapping(publish_job["permissions"]) == {"id-token": "write"}
        assert _mapping(publish_job["environment"])["name"] == "pypi"
        publish_steps = _steps(publish_job)
        assert all("run" not in step for step in publish_steps)
        assert [step.get("uses") for step in publish_steps] == [
            DOWNLOAD_ARTIFACT_ACTION,
            PYPI_ACTION,
        ]
        download_with = _mapping(publish_steps[0]["with"])
        assert download_with["path"] == "dist"
        assert download_with["name"] == artifact_with["name"]
        assert _mapping(publish_steps[1]["with"])["packages-dir"] == "dist"
    assert "  release:" in publish
    assert "types: [published]" in publish
    assert "\n  push:" not in publish
    assert "\n  workflow_dispatch:" not in publish
    assert publish.count("id-token: write") == 2
    assert publish.count(PYPI_ACTION) == 2
    assert publish.count("name: pypi") == 2
    assert "--kind final" in publish
    assert "--kind rc" in publish
    assert "github.event.release.prerelease" in publish
    assert publish.count("github.event.release.tag_name") >= 4
    assert publish.count("--expected-commit") == 2
    assert publish.count("--dist dist") == 4
    assert publish.count("python -m scripts.rebuild_wheel_from_sdist") == 2
    assert publish.count("build/sdist-wheel/*.whl") == 2
    assert publish.count(UPLOAD_ARTIFACT_ACTION) == 2
    assert publish.count(DOWNLOAD_ARTIFACT_ACTION) == 2


def test_continuous_and_build_only_workflows_cannot_upload() -> None:
    for name in ("py-ci.yml", "nightly.yml", "release-build.yml"):
        workflow = _workflow(name)
        document = _workflow_document(name)
        jobs = _mapping(document["jobs"])
        assert set(jobs) == {"py" if name == "py-ci.yml" else "build"}
        steps = _steps(next(iter(jobs.values())))
        assert any(step.get("name") == "Verify distribution identity" for step in steps)
        assert any(step.get("name") == "Rebuild wheel from sdist" for step in steps)
        assert any(
            "build/sdist-wheel/*.whl" in str(step.get("run", "")) for step in steps
        )
        assert "id-token: write" not in workflow
        assert "gh-action-pypi-publish" not in workflow

    nightly_events = _mapping(_workflow_document("nightly.yml")["on"])
    assert set(nightly_events) == {"schedule", "workflow_dispatch"}
    assert isinstance(nightly_events["schedule"], list)
    nightly = _workflow("nightly.yml")
    assert "retention-days: 14" in nightly
    assert 'id: source\n        run: echo "sha=$(git rev-parse HEAD)"' in nightly
    assert "bioextract-nightly-${{ steps.source.outputs.sha }}" in nightly
    assert "bioextract-nightly-${{ github.sha }}" not in nightly
    release_build = _workflow("release-build.yml")
    release_events = _mapping(_workflow_document("release-build.yml")["on"])
    assert set(release_events) == {"push", "workflow_dispatch"}
    assert "--expected-commit" in release_build
    assert "--dist dist" in release_build


def test_ci_installs_project_and_verifies_exact_wheel() -> None:
    ci = _workflow("py-ci.yml")
    ci_events = _mapping(_workflow_document("py-ci.yml")["on"])
    assert set(ci_events) == {"push", "pull_request"}
    assert "pdm sync -G dev\n" in ci
    assert "--no-self" not in ci
    assert "scripts/verify_installed_wheel.py" in ci
    assert "python -m scripts.rebuild_wheel_from_sdist" in ci
    assert "--force-reinstall --no-deps build/sdist-wheel/*.whl" in ci
