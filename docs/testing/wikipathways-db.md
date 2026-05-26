# WikiPathwaysDb Test Plan

## Scope

The first WikiPathwaysDb tests verify GMT parsing, metadata extraction,
species filtering, enrichment input extraction, single and grouped selections,
tidy writing, and error handling.

The tests do not cover online WikiPathways APIs, GPML parsing, pathway topology,
or enrichment p-value calculations.

## Fixtures

Use a tiny GMT file that preserves the real WikiPathways shape:

```text
Glutathione metabolism%WikiPathways_20260510%WP100%Homo sapiens	https://www.wikipathways.org/instance/WP100	2687	2678
Alanine and aspartate metabolism%WikiPathways_20260510%WP106%Homo sapiens	https://www.wikipathways.org/instance/WP106	2806	435
Mouse pathway%WikiPathways_20260510%WP1%Mus musculus	https://www.wikipathways.org/instance/WP1	123
```

## Unit Tests

- `from_gmt()` accepts an existing GMT file.
- `from_gmt()` rejects missing and oversized files.
- Malformed pathway headers raise a targeted `ValueError`.
- `extract_pathway()` parses pathway name, collection, version, identifier,
  species, URL, and gene count.
- `extract_term2gene()` emits `WikiPathwaysId, GeneId`.
- `extract_term2name()` emits pathway display metadata.
- `species=` filters pathway metadata and filters `term2gene` through retained
  pathway IDs.
- `select_ids()` trims input IDs and drops blanks.
- `extract_mapping()` returns selected Entrez IDs joined to pathway metadata.
- `extract_unmapped_input_ids()` reports IDs not present in the GMT gene sets.
- `select_groups()` preserves `GroupId`.
- `build_tidy().write()` writes flat parquet files and optional manifest.
- Resource count limits are enforced.

## Real-Data Smoke

Use:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/wikipathways/gmt/2026-05-10/raw/Homo_sapiens/wikipathways-20260510-gmt-Homo_sapiens.gmt
```

Check that:

- `extract_pathway()` is non-empty.
- `extract_term2gene()` is non-empty.
- `extract_term2name()` is non-empty.
- tidy writing succeeds under `/tmp`.
