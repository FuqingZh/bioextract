# UniprotDb Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

The UniProtDb test suite verifies lightweight construction, raw
`idmapping_selected` parsing, taxid filtering, single parquet writing, hive
dataset reading, single parquet reading, eggNOG xref flat-file extraction, and
Swiss-Prot subcellular location flat-file extraction, and schema/error
handling.

It does not test online UniProt services or full 9 GB all-taxid exports.

## Fixtures

Use tiny raw `.tab` and `.tab.gz` fixtures with all 22 selected idmapping
columns and mixed `TaxId` values:

```text
P04637	P53_HUMAN	7157	...	9606	...
Q9Y243	AKT3_HUMAN	10000	...	9606	...
P31750	AKT1_MOUSE	11651	...	10090	...
```

## Unit Tests

- `from_files()` accepts raw `.tab.gz`.
- `from_files()` accepts raw `.tab`.
- `from_files()` accepts single `mapping.parquet`.
- `from_files()` accepts hive parquet dataset directories.
- missing, unsupported, empty hive, and oversized inputs are rejected.
- `with_taxids()` normalizes taxid values and rejects empty values.
- `extract_mapping()` filters by selected taxids.
- `write_tidy()` emits canonical `mapping.parquet`.
- `level_compression` is passed through to zstd parquet writing.
- all-taxid write requires `should_allow_all=True`.
- non-empty output directories follow `policy_existing`.
- `validate_schema()` reports missing required columns.
- `.dat(.gz)` snapshots can emit eggNOG xref tidy output.
- `.dat(.gz)` snapshots can extract subcellular location comments.
- subcellular fixtures cover simple locations, multiple locations, ECO
  evidence, no-evidence comments, and `Note=` text.
- `write_subcellular_location_tidy()` emits canonical `data.parquet`.

## Real-Data Smoke

For the local 2026_01 snapshot:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/uniprot/idmapping/2026_01/raw/idmapping_selected.tab.gz
```

Keep smoke tests conservative because the file is about 9.4 GB compressed:

- `from_files()` should succeed.
- optional `validate_schema()` can be run as a preflight.
- taxid-scoped extraction can be tested on a known small taxid if needed.

For Swiss-Prot knowledge-base validation, the local flat-file snapshot is:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/uniprot/kb/2026_01/raw/knowledgebase/complete/uniprot_sprot.dat.gz
```

The eggNOG xref tidy path should verify:

- `.dat.gz` parsing succeeds
- `mapping.parquet` is written
- `EggnogOgId` and `EggnogLevel` are populated from `DR   eggNOG;` records

The subcellular location tidy path should verify:

- `.dat.gz` parsing succeeds
- `data.parquet` is written
- distinct primary accessions and row counts are stable for the snapshot
- evidence coverage is non-zero
- location text and `Note=` text are populated in a deterministic preview
