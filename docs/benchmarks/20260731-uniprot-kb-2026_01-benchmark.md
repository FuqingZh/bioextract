# UniProtKB 2026_01 Publication Baseline

Date: 2026-07-31
Status: legacy pre-convergence metadata-v3 DuckDB baseline

This artifact predates the supported metadata-v1 contract and is retained only
as historical measurement evidence. Current metadata-v1-only readers reject
it; it is not an accepted current publication. A formal metadata-v1 rebuild
benchmark remains pending.

## Artifact

The reviewed UniProtKB/Swiss-Prot 2026_01 snapshot was published to:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/uniprot/kb/2026_01/tidy/uniprot.duckdb
```

The artifact was built from `bioextract` main commit
`038c8215a1ca2f1daae58514de62a8fb06067395`. It is 1,498,165,248 bytes and
has SHA-256
`98ef664b1426af609b09bf0c976a5086bd31432c9f38ae12989f3a0dd1d3cabe`.

Publication completed in 10:39.75 wall-clock time, used 746.27 seconds of user
CPU and 65.52 seconds of system CPU at 126% CPU, and peaked at 12,280,916 KiB
resident memory. No swaps occurred.

These measurements are an observational baseline for this snapshot and host,
not a performance service-level agreement.

## Inputs And Metadata

The three caller-declared source roles were:

| Role | Bytes | SHA-256 |
| --- | ---: | --- |
| UniProtKB DAT | 692,563,345 | `bb3815e7b6445566ad9c8479f659033aa2115ed3cf2b06e61ae37c1dabc60438` |
| Canonical FASTA | 93,457,057 | `5ba5cb332fc7794ab1c02075a79c8b3d95b573f9b244a38bb53558172e1f9b7b` |
| Varsplic FASTA | 8,575,451 | `ce8bb549175f6901722b0e45a7d8685edffd4d36c49ceaf9bffb1fbc075044d4` |

Embedded metadata records:

- metadata schema version `3`;
- resource `uniprot`;
- resource schema `uniprot-knowledgebase-duckdb-v1`;
- source profile `uniprotkb-flat-file-v1`;
- release `2026_01`, with release source `caller`;
- molecular-weight validation model `legacy-expasy`; and
- `validation_issue_count` 0.

The `_bioextract.source_file` inventory contains all three source roles with
the sizes shown above.

## Publication Inventory

The publication contains 16 biological tables and 38,259,122 rows:

| Relation | Rows | Relation | Rows |
| --- | ---: | --- | ---: |
| `protein` | 574,627 | `protein_accession` | 823,596 |
| `protein_sequence` | 615,960 | `protein_isoform` | 70,077 |
| `protein_isoform_identifier` | 70,446 | `protein_sequence_variation` | 53,476 |
| `protein_isoform_variation` | 62,535 | `protein_name` | 1,192,108 |
| `gene_name` | 1,132,627 | `protein_ec_number` | 305,721 |
| `protein_go_annotation` | 3,358,100 | `protein_cross_reference` | 18,465,377 |
| `protein_comment` | 2,802,336 | `protein_subcellular_location` | 474,270 |
| `protein_keyword` | 3,839,520 | `protein_identifier` | 4,418,346 |

## Domain Acceptance

Validated reopening with `UniProtDatabase.from_duckdb()` and native read-only
SQL both passed. Write SQL through the native connection was rejected because
the artifact is read-only.

Formal biological probes also passed:

- secondary accession `P18556` resolves to both `P68744` and `P68745`;
- protein `Q6T412` has sequence length 421 and CRC64
  `EE7D1FA88E010B94`; and
- protein `P04637` has nine isoform products.

The completed publication left no staging file, DuckDB WAL, relation spool, or
validation-index residue.

## Warm Query Observations

Warm queries against the CephFS artifact produced:

| Operation | Runs | Median | Approximate p95 |
| --- | ---: | ---: | ---: |
| DuckDB point lookup | 100 | 9.765 ms | 12.277 ms |
| DuckDB GO join | 50 | 18.051 ms | 21.097 ms |
| Domain isoform query | 20 | 51.244 ms | 58.960 ms |

These warm-query measurements are observations from this snapshot and host,
not a performance service-level agreement.
