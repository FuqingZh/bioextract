# KEGG Metabolic 2026-07 Publication Baseline

Date: 2026-07-30
Status: legacy pre-convergence DuckDB baseline

This artifact predates the supported metadata-v1 contract and is retained only
as historical measurement evidence. It is not an accepted current publication;
a formal metadata-v1 rebuild benchmark remains pending.

## Artifact

The complete local KEGG metabolic snapshot was published to:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/kegg/metabolic/2026-07/tidy/kegg.duckdb
```

The source contained 4,112 files totaling 467,360,027 bytes. Publication with
source SHA-256 calculation completed in 2:08.57 wall-clock time, peaked at
1,476,904 KiB resident memory, and produced a 17,313,792-byte DuckDB.

These measurements are an observational baseline for this snapshot and host,
not a performance service-level agreement.

## Publication Inventory

The publication contains 25 biological tables and 444,929 rows:

| Relation | Rows | Relation | Rows |
| --- | ---: | --- | ---: |
| `compound` | 19,619 | `compound_cross_reference` | 71,540 |
| `compound_name` | 33,352 | `compound_pubchem` | 19,443 |
| `compound_reaction` | 52,454 | `enzyme` | 8,343 |
| `enzyme_cross_reference` | 31,941 | `enzyme_ko` | 10,084 |
| `enzyme_name` | 31,977 | `enzyme_replacement` | 1,441 |
| `module` | 573 | `module_compound` | 3,538 |
| `module_definition_node` | 5,882 | `module_ko` | 4,171 |
| `module_pathway` | 1,681 | `module_reaction_step` | 3,187 |
| `reaction` | 12,459 | `reaction_class` | 16,626 |
| `reaction_cross_reference` | 6,846 | `reaction_enzyme` | 10,752 |
| `reaction_ko` | 12,276 | `reaction_module` | 3,182 |
| `reaction_name` | 9,572 | `reaction_participant` | 54,162 |
| `reaction_pathway` | 19,828 |  |  |

All 4,112 source records retain their original display paths and 64-character
SHA-256 values. Metadata records release `2026-07`, the eleven discovered
source-role capabilities, and `passed_with_warnings`.

The two warnings are expected source-integrity observations: reactions
`R02417` and `R13604` reference missing compound `C23109`. Both are persisted
as `foreign_key_violation` records in `_bioextract.validation_issue`.

## Domain Acceptance

Validated reopening with `KEGGDatabase.from_duckdb()` and native read-only SQL
both passed. Write SQL through `connect()` was rejected by DuckDB.

All nine public namespaces resolved against the formal artifact:

| Namespace | Probe | Matched anchors | Selected reactions |
| --- | --- | ---: | ---: |
| `kegg_compound` | `C00001` | 1 | 4,378 |
| `chebi` | `CHEBI:15377` | 1 | 4,378 |
| `pubchem` | `3303` | 1 | 4,378 |
| `kegg_reaction` | `R00001` | 1 | 1 |
| `rhea` | `RHEA:22455` | 1 | 1 |
| `ec` | `3.6.1.10` | 1 | 1 |
| `ko` | `K01457` | 1 | 1 |
| `kegg_module` | `M00532` | 1 | 12 |
| `kegg_pathway` | `map00220` | 33 | 33 |

Grouped selection retained group isolation: a known ChEBI ID matched only its
own group, while `CHEBI:999999999` produced `not_found` only in the missing-ID
group. Supplying all 4,171 published KOs evaluated all 573 modules as complete,
including recursive module references.

Shared identifiers joined by direct string equality without casts or prefix
construction:

- all 17,192 KEGG ChEBI cross-reference rows joined to
  `ChEBI compound.chebi_id`;
- all 6,846 KEGG Rhea cross-reference rows joined to
  `Rhea reaction.accession`.

## Warm Query Observations

Warm in-process queries against the CephFS artifact produced:

| Operation | Runs | Median | Approximate p95 |
| --- | ---: | ---: | ---: |
| DuckDB compound point lookup | 100 | 1.77 ms | 3.28 ms |
| DuckDB compound-to-pathway join | 100 | 7.94 ms | 11.06 ms |
| Domain `select_ids()` plus pathway extraction | 20 | 66.08 ms | 128.58 ms |

The domain measurement includes selection construction, validation, lineage,
and extraction overhead. It is not directly comparable to a single raw SQL
statement.
