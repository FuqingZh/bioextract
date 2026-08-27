# EggNOGDatabase Test Standard

Version: v1.1
Date: 2026-08-27
Status: current

## Scope

The EggNOGDatabase test standard covers:

- lightweight construction
- COG function lookup parsing
- single and grouped eggNOG protein selection
- unmapped reporting
- gzip warning and temporary-storage cleanup

It does not cover:

- emapper runtime behavior
- DIAMOND alignment
- KEGG/module/domain annotation parity
- enrichment statistics

## Unit Tests

- plain SQLite is queried directly without a warning.
- `from_sqlite()` detects gzip transport by content and warns exactly once at
  the caller location.
- `EggnogSelection.mappings()` expands OGs and COG categories correctly.
- `select_ids()` returns filtered rows in the fixed `eggnog_protein` namespace.
- `select_groups()` resolves one globally unique ID set and preserves `group_id`.
- unmapped IDs are reported correctly.
- pipe-bearing caller text remains a trimmed exact eggNOG protein input; it is
  not parsed as a UniProt representation in single or grouped selection.
- mapping and unmatched extractors reuse one Selection cache.
- gzip SQLite snapshots use `temp_dir` only as scratch storage, clean up after
  success and failure, and create no persistent decompressed copy.
- the database-level eager extractor and publication writers are absent.
- invalid `namespace` raises targeted `ValueError`.

## Real-Data Validation

Real-snapshot validation must verify the exact schema, a stable row count for
the selected snapshot, and cleanup of the temporary SQLite workspace. Keep
this outside the default pytest suite because the compressed database and
decompressed workspace are large.

The observed eggNOG 5.0.2 sizes and resource usage are recorded in the
[eggNOG 5.0.2 benchmark](../benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md).
