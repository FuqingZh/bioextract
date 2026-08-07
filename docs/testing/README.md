# bioextract Test Standard

Version: v1.0
Date: 2026-07-31
Status: current

## Purpose

bioextract tests protect four different boundaries. Keep those boundaries
explicit so a failing test identifies the owner of the regression and focused
development does not require running every publication pipeline.

## Test Layers

### Unit

`tests/unit/` verifies one component in isolation: parsers, identifier
normalization, equation or sequence algorithms, and shared frame construction.
A unit test may use an in-memory frame or a small temporary file when that is
the component's native input. It does not build and reopen a complete
publication.

### Contract

`tests/contract/` verifies stable caller and artifact contracts:

- top-level and resource-package imports;
- method signatures, documented behavior, and stable exceptions;
- DuckDB and Parquet schemas;
- embedded metadata and provenance;
- capability and validation semantics;
- atomic replacement and corrupt-publication rejection.

Contract tests may create a minimal temporary DuckDB or Parquet file. File I/O
alone does not make a test an integration test.

### Integration

`tests/integration/` verifies that real collaborating components work together
using compact, local, official-format fixtures:

```text
source fixture
    -> resource parser
    -> publication build
    -> reopen
    -> select and extract
```

Integration tests use the real Polars, DuckDB, SQLite, archive, and filesystem
paths. They remain hermetic, deterministic, and part of the canonical gate.

### External Snapshot Smoke

`tests/smoke/` is reserved for read-only checks against an installed package or
host-owned official publications such as CephFS snapshots. Smoke tests require
explicit invocation and never run as part of the default test suite.

A test using only `tmp_path` and generated local data is not a smoke test even
when its name says `smoke`; it belongs to unit, contract, or integration
according to the boundary it crosses.

## Classification Rules

Classify by the largest boundary exercised, not by implementation detail:

1. complete build/reopen/query flow: integration;
2. stable public or publication invariant: contract;
3. one parser, normalization, or algorithm component: unit;
4. installed artifact or external snapshot: smoke.

The directory is the source of truth. Do not add `unit`, `contract`, or
`integration` markers that duplicate it. Markers are reserved for orthogonal
requirements such as an external snapshot.

## Fixtures And Support

- Prefer explicit fixtures over autouse setup.
- Keep resource fixture factories in the nearest resource test directory.
- Put stable official-format samples under `tests/fixtures/<resource>/`.
- Keep malformed-input builders beside the tests that explain the malformed
  contract.
- Put only format-neutral, genuinely shared mechanics in `tests/_support/`.
- Keep root `conftest.py` thin; do not turn it into a repository-wide fixture
  registry.
- Use `tmp_path`; never write tests to the real resource tree or CephFS.

## Commands

```console
pdm run test-unit
pdm run test-contract
pdm run test-integration
pdm run test
pdm run check
```

`pdm run test` and `pdm run check` cover all hermetic layers. Resource-specific
development may pass a narrower path to pytest.

External smoke is explicit:

```console
pdm run test-smoke
```

Configure one or more DuckDB publications with the platform path separator:

```console
BIOEXTRACT_SMOKE_DUCKDBS=/path/a.duckdb:/path/b.duckdb \
  pdm run test-smoke
```

The smoke command skips when no host publication is configured; the hermetic
suite must not skip because of external state.

## Resource Limits

Repository-owned validation defaults DuckDB, Polars, Rayon, BLAS, and OpenMP
worker pools to a common ceiling of four threads. The launcher applies that
environment before each `lock-check`, format check, lint, typecheck, or pytest
child starts. External-publication smoke defaults to two. Override the common
ceiling when a smaller shared-host budget is required:

```console
BIOEXTRACT_TEST_THREADS=1 pdm run check
```

`BIOEXTRACT_TEST_THREADS` accepts any positive integer. Invalid budgets fail
before the validation child starts. An unset native-pool variable receives the
budget, a smaller positive caller value is preserved, and a larger positive
value is clamped. Invalid native-pool values also fail before validation.

The ceiling bounds known native worker pools; helper threads mean it is not an
exact cap on total process threads or NLWP. It belongs to repository validation
only and does not change bioextract runtime defaults for real publication
builds. Host-wide validation admission and process-tree containment remain the
responsibility of the host runner.

## Naming

Name tests for the behavior and expected outcome. Use `round_trip` for a
minimal local source-to-query flow. Reserve `smoke` for installed-package or
external-publication checks.

Avoid encoding the implementation module name when the stable behavior is a
better identifier.

## Gate Contract

- Focused layer tests should pass before the complete suite.
- The complete suite is required for public API, publication schema, or
  cross-resource changes.
- A structural migration must preserve the pre-migration collection count
  unless a mixed test is deliberately split and the increase is documented.
- Full release snapshots and performance benchmarks remain outside default
  pytest and follow their resource testing or benchmark document.
