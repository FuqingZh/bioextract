# bioextract Agent Guide

## Project Boundary

- `bioextract` turns local biological-database snapshots into stable,
  provenance-aware extraction outputs. It does not download source releases.
- Use `docs/README.md` as the authority map before changing a public API,
  identifier rule, output schema, or publication behavior.
- Keep resource-specific parsing and biological semantics in the owning package
  under `src/bioextract/`.
- Real snapshot builds and external storage writes require explicit task scope;
  tests should use temporary local fixtures.
- Before delivery, run the non-mutating code gate used by CI:

  ```console
  pdm run check
  ```

- `pdm run format` and `pdm run lint-fix` mutate files; use them only on the
  intended change set.

## AO Delivery

This repository has opted into the accepted user-level AO service as
`bioextract`. For conversation-authorized implementation intended to cross a
pull-request boundary, verify AO health and start a task-specific worker before
creating the implementation branch or pull request. If a pull request already
exists, mark it ready for review if it is a draft, then restore its owning
worker or claim it with `--no-takeover`. Ready-for-review is only an AO claim
prerequisite; leave merge and risk decisions to the user. If AO is unavailable,
use an isolated worktree and report that fallback.

Classify AO observations by state owner before diagnosing them:

- sandbox state: paths and processes visible only inside the current agent
  sandbox;
- worker state: the AO-created worktree, pane, environment, and Codex home;
- daemon state: the persistent AO service, database, and project readback; and
- host state: the user service manager, filesystem, credentials, and installed
  binaries outside the worker boundary.

A mismatch between these states is diagnostic evidence, not proof that the
host is broken. Verify the state through its owning context and use the AO
runbook's repair procedure before changing persistent host configuration.

Use the accepted AO diagnosis states exactly:

- a failure observed only in the sandbox is `indeterminate`;
- active host service plus AO `ready`/`running` readback and a passing
  `healthz` probe is `daemon ready`;
- repeated failure from the authoritative host context is `unavailable`; and
- an AO doctor external integration or authentication failure is `delivery
  degraded`, not daemon unavailability. Core doctor failures remain evidence
  about daemon or host readiness and must be diagnosed by their owning state.
