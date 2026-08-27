from __future__ import annotations

import re
from pathlib import Path
from typing import cast

import pytest
import yaml

REPOSITORY_ROOT = Path(__file__).parents[3]
WORKFLOW_ROOT = REPOSITORY_ROOT / ".github" / "workflows"
PYPI_ACTION = "pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33"


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
    assert set(jobs) == {"final", "rc"}
    for job_name, prerelease_expression in (
        ("final", "!github.event.release.prerelease"),
        ("rc", "github.event.release.prerelease"),
    ):
        job = _mapping(jobs[job_name])
        assert prerelease_expression in cast("str", job["if"])
        assert _mapping(job["permissions"]) == {
            "contents": "read",
            "id-token": "write",
        }
        assert _mapping(job["environment"])["name"] == "pypi"
        steps = _steps(job)
        upload_steps = [step for step in steps if step.get("uses") == PYPI_ACTION]
        assert len(upload_steps) == 1
        verification_steps = [
            step for step in steps if str(step.get("name", "")).startswith("Verify ")
        ]
        assert any(
            "--dist dist" in str(step.get("run", "")) for step in verification_steps
        )
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
    assert publish.count("--dist dist") == 2


def test_continuous_and_build_only_workflows_cannot_upload() -> None:
    for name in ("py-ci.yml", "nightly.yml", "release-build.yml"):
        workflow = _workflow(name)
        document = _workflow_document(name)
        jobs = _mapping(document["jobs"])
        assert set(jobs) == {"py" if name == "py-ci.yml" else "build"}
        steps = _steps(next(iter(jobs.values())))
        assert any(step.get("name") == "Verify distribution identity" for step in steps)
        assert "id-token: write" not in workflow
        assert "gh-action-pypi-publish" not in workflow

    nightly_events = _mapping(_workflow_document("nightly.yml")["on"])
    assert set(nightly_events) == {"schedule", "workflow_dispatch"}
    assert isinstance(nightly_events["schedule"], list)
    assert "retention-days: 14" in _workflow("nightly.yml")
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
