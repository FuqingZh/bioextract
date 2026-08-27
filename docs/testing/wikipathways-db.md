# WikiPathwaysDatabase Test Standard

Version: v1.1
Date: 2026-08-27
Status: current

## Scope

The WikiPathwaysDatabase tests verify resource-local source resolution,
whole-dataset GMT validation, metadata extraction, species filtering,
enrichment input extraction, single and grouped selections, tidy writing, and
error handling.

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

- `from_gmt()` accepts a literal scalar, a literal sequence, and glob
  expressions.
- `**` glob expressions resolve files in nested directories recursively.
- Resolution freezes normalized actual files in deterministic order.
- `glob=False` treats scalar and sequence entries literally.
- Empty input, unmatched patterns, missing paths, directories, and non-files
  are rejected.
- Repeated literals, overlapping globs, symlink aliases, and hard links that
  reach the same physical file are rejected.
- The complete resolved set requires one official Collection and exposes the
  one Version derived from that Collection.
- Collection values without the `WikiPathways_` prefix or its non-empty Version
  suffix are rejected before release provenance is published.
- Each resolved GMT must contain a non-empty pathway record; both a single
  empty file and an empty file alongside a valid file raise a path-specific
  error.
- Duplicate `WikiPathwaysId` values within one file or across files are
  rejected.
- Malformed pathway headers raise a targeted `ValueError`.
- `pathways()` parses pathway name, collection, version, identifier,
  species, URL, and gene count.
- `pathway_genes()` emits `wiki_pathways_id, gene_id`.
- `pathway_names()` emits pathway display metadata.
- A single GMT may contain multiple species.
- `species=` filters pathway metadata and `term2gene` by row content after
  whole-dataset validation.
- `select_ids()` trims input IDs and drops blanks.
- pipe-bearing caller text remains a trimmed exact Entrez lookup/unmatched
  value in source and reopened single/grouped selection.
- `mappings()` returns selected Entrez IDs joined to pathway metadata.
- `unmatched_ids()` reports IDs not present in the GMT gene sets.
- `select_groups()` preserves `group_id`.
- `write_duckdb()` writes pathway and pathway-gene relations in one file.
- `_bioextract` contains every resolved actual source under deterministic
  unique logical names, including files with no retained rows, plus row-count
  provenance.
- `_bioextract.metadata` records the common content-derived Version as
  `bioextract.release_version` with
  `bioextract.release_version_source=official_metadata`; the publication test
  asserts both values.

## DuckDB Reopen

- Build a representative species-scoped GMT source, publish it, reopen it with
  `from_duckdb()`, and compare `pathways()`, `pathway_genes()`,
  `pathway_names()`, single-selection, and grouped-selection outputs with the
  source-backed handle.
- `connect()` returns distinct caller-owned native read-only connections,
  permits arbitrary read SQL, and rejects writes. Source-backed `connect()`
  raises `CapabilityError` without changing existing source extraction errors.
- Reject metadata schema v1 and require metadata schema v2, source profile
  `wikipathways-gmt-v1`, and the current resource schema. Reject incompatible
  resource identity, release identity, exact physical table/view inventory,
  table roles, column provenance, and physical column schemas.
- Reopen validation accepts non-negative recorded biological row counts without
  recounting or scanning the biological tables.
- Reopened handles and cached selection terminals reject vanished or atomically
  replaced files. Reopened tidy datasets and handles cannot be republished.

## Real-Data Smoke

Use:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/wikipathways/gmt/2026-05-10/raw/Homo_sapiens/wikipathways-20260510-gmt-Homo_sapiens.gmt
```

Check that:

- `pathways()` is non-empty.
- `pathway_genes()` is non-empty.
- `pathway_names()` is non-empty.
- tidy writing succeeds under `/tmp`.
