# UniprotDb Test Plan

## Scope

The first UniProtDb test suite verifies lightweight construction, raw
`idmapping_selected` parsing, taxid filtering, hive parquet writing, hive
dataset reading, single parquet reading, and schema/error handling.

It should not test online UniProt services, `.dat` parsing, or full 9 GB
all-taxid exports.

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
- `write_tidy()` defaults to hive partitioning.
- hive output can be read back through `from_files()`.
- all-taxid write requires `should_allow_all=True`.
- non-hive write is allowed only for one selected taxid.
- `validate_schema()` reports missing required columns.

## Real-Data Smoke

For the local 2026_01 snapshot:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/uniprot/idmapping/2026_01/raw/idmapping_selected.tab.gz
```

Keep smoke tests conservative because the file is about 9.4 GB compressed:

- `from_files()` should succeed.
- optional `validate_schema()` can be run as a preflight.
- taxid-scoped extraction can be tested on a known small taxid if needed.
