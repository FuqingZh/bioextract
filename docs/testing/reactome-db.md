# ReactomeDatabase Test Standard

Version: v1.1
Date: 2026-08-17
Status: current

## Scope

The ReactomeDatabase suite verifies local official-format parsing, explicit
lowest/all-level capability semantics, selection behavior, species scoping,
enrichment projections, hierarchy validation, metadata-v2 publication, and
DuckDB reopen parity. It does not test online Reactome services or statistical
enrichment.

## Fixtures

Keep compact fixtures headerless and tab-separated. Mapping fixtures use the
six columns:

```text
uniprot_id
reactome_pathway_id
reactome_url
pathway_name
evidence_code
species
```

Provide separate `UniProt2Reactome.txt` and
`UniProt2Reactome_All_Levels.txt` fixtures. The all-level fixture must not be
passed through the lowest-level role. Include:

- an ancestor mapping present only at all levels;
- an exact duplicate row;
- two rows differing only by `evidence_code`;
- an embedded double quote in a reaction-shaped six-field fixture; and
- a missing pathway metadata ID for publication-warning coverage.

Tests use `tmp_path` and never write to the formal resource tree or CephFS.

## Unit Tests

- The shared mapping-family reader disables CSV quoting and preserves literal
  quotes.
- Short, extra-field, null, and empty-field records fail closed.
- Exact duplicate six-column rows are removed; evidence-distinct rows remain.
- Pathway and relation readers preserve their declared schemas.

## Integration Tests

### Construction and query

- `from_files()` accepts any non-empty combination of the four P1 roles.
- `uniprot_all_levels` is independent from `uniprot_mapping`.
- `release_version` trims caller input and rejects an empty value.
- Missing paths and non-files fail at construction.
- Default `pathway_mappings()`, `pathway_genes()`, `select_ids()`, and
  `select_groups()` are lowest-level UniProt behavior.
- `pathway_level="all_levels"` is explicit for whole-resource and selected
  queries; absent capability raises `CapabilityError` or the source-handle
  `ValueError` at the feature boundary.
- Invalid namespace, target, and level values raise `ValueError`.
- Evidence-distinct selected rows survive; `pathway_genes()` deduplicates only
  pathway/UniProt pairs.
- Grouped mapping and unmatched output retain their existing shapes.
- Species filtering occurs before enrichment deduplication.

### Publication and reopen

For lowest-only, all-only, pathway-only, relation-only, and combined fixtures:

- `build_tidy().frames` contains only the available canonical role names:

  ```text
  uniprot_pathway_lowest_level
  uniprot_pathway_all_level
  pathway
  pathway_relation
  ```

- `write_duckdb()` publishes exactly those biological tables plus the five
  metadata-v2 tables.
- Metadata identifies `reactome-mapping-files-v2` and
  `reactome-mapping-v0.2`.
- `_bioextract.source_file` contains the exact role inventory.
- The old `protein_pathway` relation is absent.
- `release_version` survives publication, inspection, and reopen.
- Missing pathway metadata creates a visible
  `missing_pathway_metadata` warning without dropping mapping rows.
- When all three UniProt pathway roles are present, all-level keys equal the
  reflexive lowest-level hierarchy closure; a mismatch fails publication.
- Reopen derives capabilities from the unique source-role inventory and rejects
  forged roles, wrong profile/version, extra or missing physical tables,
  physical schema changes, and v0.1 publications.
- Reopen trusts non-negative recorded biological row counts and does not recount
  them during validation.
- Source and reopened whole-resource, species-scoped, selected, grouped,
  unmatched, enrichment, and relation results agree.
- Reopened `connect()` calls are distinct read-only caller-owned connections,
  and atomic replacement invalidates the old handle.

## Bounded v96 Smoke

External snapshot smoke is opt-in and read-only. Use the concrete raw subtree
only after resolving it explicitly; do not recursively scan `/cephfs_data`.
The smoke writes any temporary publication under a unique `/tmp` directory.

The v96 checks record, at minimum:

```text
lowest-level rows                    322435
all-level rows                       934947
lowest-level unique pathway pairs   317978
all-level unique pathway pairs      917206
human lowest-level rows              53996
human all-level rows                 159114
human unique UniProt IDs              12136 at both levels
closure derived-only keys                 0
closure official-only keys                 0
```

The smoke also proves that the 13 all-level rows for
`R-SCE-9865878` remain queryable and that the corresponding publication
warning is visible. It records runtime and output hashes only when the bounded
run actually completes; fixture results are not substituted for skipped
external evidence.

## Commands

```console
BIOEXTRACT_TEST_THREADS=1 pdm run pytest \
  tests/unit/reactome \
  tests/integration/reactome \
  tests/contract/resources/reactome \
  -q

BIOEXTRACT_TEST_THREADS=1 pdm run check
```

The complete gate is required for this public API and publication-schema
change. Formal CephFS replacement, release catalog admission, and package
release are outside this test standard.
