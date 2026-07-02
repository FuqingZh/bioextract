# InterProDb Test Plan

## Scope

The InterProDb test plan covers:

- lightweight construction
- `protein2ipr.dat.gz` parsing
- optional XML enrichment
- full mapping extraction
- single and grouped UniProt selection
- unmapped reporting
- lazy tidy writing of one canonical parquet

It does not cover:

- InterProScan runtime behavior
- sequence search
- external service calls
- enrichment statistics

## Unit Tests

- `from_mapping_files()` accepts required and optional files.
- `extract_mapping()` returns the normalized mapping columns.
- XML enrichment fills `InterProType` and `MemberDb` when present.
- `select_ids(..., kind_input_id="uniprot")` returns filtered rows.
- `select_groups(..., kind_input_id="uniprot")` preserves `GroupId`.
- unmapped IDs are reported correctly.
- `write_tidy()` writes canonical `mapping.parquet`.
- invalid `kind_input_id` raises targeted `ValueError`.

## Real-Data Validation

Validated inputs:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/interpro/mapping/108.0/raw/protein2ipr.dat.gz
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/interpro/mapping/108.0/raw/interpro.xml.gz
```

Validated output:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/interpro/mapping/108.0/tidy/mapping.parquet
```

Observed publication result on 2026-06-08:

- parquet size: about 13 GB
- row count: 1,175,529,272
- schema:
  - `UniProtId`
  - `InterProId`
  - `InterProName`
  - `InterProType`
  - `MemberDb`
  - `MemberDbId`
  - `Start`
  - `End`

Observed runtime characteristics:

- `POLARS_MAX_THREADS=8`
- thread count stayed around 31
- RSS was observed up to about 17 GB
- the invalid earlier partial parquet was removed before rerun

Smoke validation after write:

- parquet schema could be scanned
- row count could be read via `pl.scan_parquet(...).select(pl.len())`
- `head()` sampling returned valid rows
