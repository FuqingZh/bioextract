# EggNOGDatabase Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

The EggNOGDatabase test standard covers:

- lightweight construction
- COG function lookup parsing
- full mapping extraction
- single and grouped eggNOG protein selection
- unmapped reporting

It does not cover:

- emapper runtime behavior
- DIAMOND alignment
- KEGG/module/domain annotation parity
- enrichment statistics

## Unit Tests

- `from_sqlite()` accepts SQLite and gzip-wrapped SQLite inputs.
- `extract_mapping()` expands OGs and COG categories correctly.
- `select_ids(..., namespace="eggnog_protein")` returns filtered rows.
- `select_groups(..., namespace="eggnog_protein")` preserves `GroupId`.
- unmapped IDs are reported correctly.
- gzip SQLite snapshots are decompressed under `dir_tmp`.
- invalid `namespace` raises targeted `ValueError`.

## Real-Data Validation

Real-snapshot validation must verify the exact schema and a stable row count
for the selected snapshot. Keep this outside the default pytest suite because
the compressed database and decompressed workspace are large.

The observed eggNOG 5.0.2 sizes and resource usage are recorded in the
[eggNOG 5.0.2 benchmark](../benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md).
