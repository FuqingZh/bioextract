# Repository Delivery Boundaries

Version: v1.0
Date: 2026-08-07
Status: current

This document contains only the `bioextract`-specific additions to the active
host or orchestrator delivery policy. It is not an AO operating manual.

## Ownership

- Keep one writer for an implementation branch and worktree. A controller must
  not patch, stage, commit, or push another worker's owned worktree.
- The repository owns source code, contracts, tests, documentation, and its CI
  entrypoint. The host owns orchestrator lifecycle, session state, credentials,
  containment, repository settings, and merge policy.
- Do not persist AO commands, daemon-health rules, credentials, `CODEX_HOME`,
  or host service configuration in this repository.

## Verification

Use the smallest relevant test path while iterating. Run the complete CI gate
for a public API, publication schema, or cross-resource change:

```console
pdm run check
```

Documentation-only and narrowly scoped harness changes use checks that decide
their affected surface; they do not require the complete Python gate merely
because they are delivered through a pull request.

`BIOEXTRACT_TEST_THREADS` bounds known native-library pools used by repository
validation. It does not serialize independent workers or cap total process
threads. External snapshot smoke remains opt-in through `pdm run test-smoke`.

## Publication And Shared Storage

A code pull request does not publish to CephFS or delete a prior release. A
real snapshot build or delivery must name its exact input subtree and output
target, validate the staged artifact, close it, and replace atomically. Never
recursively scan `/cephfs_data`.
