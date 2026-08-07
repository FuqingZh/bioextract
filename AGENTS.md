# bioextract Agent Guide

## Start Here

- `docs/README.md` is the authority map. Read the linked architecture and test
  contracts before changing a public API, schema, publication boundary, or
  storage strategy.
- `bioextract` turns caller-supplied local official database snapshots into
  stable, provenance-aware domain access. It does not download releases or
  depend on `biofetch` or downstream applications.
- Keep source parsing and biological relationship semantics in the owning
  package under `src/bioextract/`. Do not add a generic query facade or move
  downstream analysis into this repository.

## Public Contract

- Database classes are the stable top-level API. Selection and result classes
  are behaviorally public but are not stable package exports.
- Python methods and parameters use `snake_case`. Generated tables, views, and
  columns use singular `snake_case`. Preserve the spelling and order of
  unchanged official two-dimensional source headers.
- Before v1.0, do not add compatibility aliases, sidecar manifests, or
  filename-derived schema identity. Embedded publication metadata is
  authoritative.

## Work And Verification

- Preserve unrelated worktree changes. Use temporary local fixtures for tests;
  real snapshots and CephFS writes require explicit task scope.
- Run focused tests while iterating. Before handing off a public-contract or
  cross-resource change, run the CI gate:

  ```console
  pdm run check
  ```

- Bound shared-host test work with `BIOEXTRACT_TEST_THREADS`. `pdm run format`
  and `pdm run lint-fix` mutate files; use them only on the intended change set.
- Keep tests in the layers documented by `docs/testing/README.md`; smoke tests
  are opt-in and must not write to the formal resource tree or CephFS.

## Safety And Delivery

- Never recursively scan `/cephfs_data`. Resolve a concrete resource subtree
  first and bound traversal, file types, file or byte size, result count, and
  concurrency.
- Publication replacement, formal CephFS delivery, and deletion of old
  artifacts are explicit release actions. Stage and replace atomically.
- For PR-bound work, read `docs/runbooks/repository-delivery.md` for this
  repository's validation, single-writer, and shared-storage boundaries.
- AO lifecycle, configuration, credentials, daemon state, merge policy, and
  `CODEX_HOME` are host-owned. Do not duplicate or persist them here.
