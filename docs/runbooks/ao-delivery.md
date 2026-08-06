# AO And Repository Delivery Norms

Version: v1.0
Date: 2026-08-04
Status: current

This document is the repository-specific companion to the calibration AO
runbook. It defines what a bioextract worker may assume and what must be read
back before a pull-request or publication handoff. Host-owned AO configuration
and credentials are not stored in this repository.

## Scope And Ownership

`bioextract` is a local-snapshot library. AO coordinates engineering work; it
does not download biological releases, publish CephFS data, or discover
catalogs. The repository remains authoritative for source code, contracts,
tests, and implementation plans. AO daemon, project, session, identity, and
worker containment state remain host-owned.

The project configuration must reference the repository `AGENTS.md` as its
worker standing instructions. It must not enable AO project `autoMerge`, add
credentials, or make a dashboard or external tracker a prerequisite for local
implementation.

## Intake And State Readback

For PR-bound implementation, read back:

```console
ao status --json
ao doctor --json
ao project get bioextract --json
ao session ls -p bioextract --json
```

When those gates hold, start a task-specific owning worker for genuinely new,
unowned PR work before creating its implementation branch or pull request.
Immediately read the new session back and hand the task to that owner through
the activity-state rules in the calibration AO runbook.

For an existing PR, mark a draft ready, then restore its owning worker. If it
has no owner, first prove every other writer quiesced, then claim it without
takeover and read the assigned owner back before sending work. Owner absence
alone does not make a ready PR unowned; ready-for-review is only a claim
prerequisite.

Classify observations by their owner: sandbox state covers only the current
agent sandbox; worker state covers the AO worktree, pane, environment, and
Codex home; daemon state covers the persistent service, database, and project
readback; host state covers the service manager, filesystem, credentials, and
installed binaries outside worker containment.

If the sandbox disagrees with these results, classify the observation by its
owner. A sandbox-only path, PID, or read-only failure is `indeterminate`; use
host-authoritative service and AO health/readiness evidence before changing
persistent host state. Active host service, AO `ready` or `running` readback,
and a passing `healthz` probe is `daemon ready`; an external integration or
authentication-only doctor failure is `delivery degraded`; repeated
authoritative core failure is `unavailable`. Only the last condition selects
the unavailable fallback. Registration and healthy configuration establish
only `runtime-ready`; do not call the repository continuation-proven until a
real anchored review or CI correction has returned to its owning worker.

When continuation is not proven, new unowned PR work uses the isolated-worktree
fallback and reports that fallback at handoff. The same reported fallback
applies when AO is authoritatively unavailable. An existing AO-owned branch,
worktree, or PR is preserved until its owner is restored or authoritative
containment and write-authority revocation are proven. The controller never
patches, stages, commits, or pushes an owner's sibling worktree.

## Worker And Branch Norms

- One worker owns one implementation branch/worktree at a time.
- Branches use `codex/<scope>` and each PR contains one independently testable
  contract or resource slice.
- The owner performs code edits, commits, pushes, CI/review observation, and
  same-scope mechanical fixes.
- Retry only bounded, idempotent transient operations; read back unknown
  external writes before retrying. Stop on head/scope changes, cancellation,
  permission or authentication failures, or retry-budget exhaustion.
- Keep dirty worktrees observable. Never force-delete a dirty worktree or claim
  descendant-process release without an empty OS-owned containment boundary.

## Verification And Handoff

The canonical repository gate is:

```console
pdm run check
```

Use focused layer commands while iterating:

```console
pdm run test-unit
pdm run test-contract
pdm run test-integration
```

Keep shared-host execution bounded with `BIOEXTRACT_TEST_THREADS=1` or `4`.
External snapshot smoke is opt-in through `pdm run test-smoke` and explicit
publication paths; it is never a reason to weaken the hermetic gate.

Before handoff, record the exact current head, focused checks, complete gate,
and any unavailable host-side check. Conversation authorization for a low-risk
implementation also authorizes native GitHub per-PR auto-merge without a second
approval, but only after exact-head required CI passes, current-head review is
clean, and no actionable review threads remain. Read all of those gates back
immediately before requesting auto-merge. An explicit user stop or a
repository-local stricter policy, or a high-risk, irreversible, permission,
security, secret, release, or compatibility decision withholds that authority
and requires escalation. This never authorizes AO project `autoMerge`.

## Publication And Shared-Storage Safety

Formal DuckDB publication builds use staging, integrity validation, connection
close, and atomic replacement. A code PR does not silently write to CephFS or
delete a prior release. Any real snapshot build or delivery names its exact
input subtree and output target; `/cephfs_data` is never recursively scanned.

## Configuration Boundary

The following are repository-owned and reviewable:

- `AGENTS.md`;
- `docs/README.md`, architecture, testing, and implementation-plans documents;
- `pyproject.toml` and `.github/workflows/` gates.

The following are host-owned and must be changed only through AO or the host
runbook:

- AO daemon, project database, session records, and worker runtime;
- credentials, `CODEX_HOME`, service units, and containment;
- GitHub permissions, repository settings, and external publication mounts.
