# bioextract Agent Guide

## Project Boundary

- `bioextract` turns local official biological-database snapshots into stable,
  provenance-aware domain access. It does not download releases or depend on
  `biofetch`, `joint_omics`, or another downstream application.
- Use `docs/README.md` as the documentation authority map. Read the shared
  domain-access and materialized-dataset contracts before changing a public
  API, identifier rule, output schema, or storage strategy.
- Keep source-specific parsing and biological relationship semantics inside
  the owning resource package under `src/bioextract/`. Do not move ordinary
  filtering or analysis into the domain API.

## Working And Verification

- Preserve unrelated worktree changes. Real snapshot builds and CephFS writes
  require explicit task scope; unit tests should use temporary local fixtures.
- Run focused tests while iterating. Before handing off a public-contract or
  cross-resource change, run the same non-mutating code gate used by CI:

  ```console
  pdm run check
  ```

- `pdm run format` and `pdm run lint-fix` mutate files; use them only on the
  intended change set.

## AO Delivery

This repository has opted into the accepted user-level AO service as
`bioextract`. For conversation-authorized implementation intended to cross a
pull-request boundary, verify `ao status --json`, `ao doctor --json`, and
`ao project get bioextract --json`, then start a task-specific worker before
creating the implementation branch or pull request.

If a pull request already exists, mark it ready for review if it is a draft,
then restore its owning worker or claim it with `--no-takeover`.
Ready-for-review is only an AO claim prerequisite.

Conversation authorization for a low-risk implementation also authorizes the
worker to request GitHub native auto-merge without a second merge
authorization, but only after required CI passes on the exact current head,
current-head review is clean, and no actionable review threads remain. Read
those gates back immediately before the request. A repository-local stricter
policy, an explicit user stop, or a high-risk, irreversible, permission,
security, secret, release, or compatibility decision withholds auto-merge and
requires escalation.

This authority applies only to GitHub's native per-pull-request auto-merge. It
does not authorize always-on AO project configuration such as `autoMerge`,
whose cancellation and state-change behavior is unproven and must remain
disabled. If AO is unavailable, use an isolated worktree and report that
fallback.

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

Do not describe this repository as continuation-proven until a real anchored
Automatic Codex Review finding has returned to its owning worker. Repository
registration and healthy configuration establish only `runtime-ready`.
