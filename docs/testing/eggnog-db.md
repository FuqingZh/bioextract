# EggnogDb Test Plan

## Scope

The EggnogDb test plan covers:

- lightweight construction
- COG function lookup parsing
- full mapping extraction
- single and grouped eggNOG protein selection
- unmapped reporting
- tidy writing through the temporary-TSV path

It does not cover:

- emapper runtime behavior
- DIAMOND alignment
- KEGG/module/domain annotation parity
- enrichment statistics

## Unit Tests

- `from_files()` accepts SQLite and gzip-wrapped SQLite inputs.
- `extract_mapping()` expands OGs and COG categories correctly.
- `select_ids(..., kind_input_id="eggnog_protein")` returns filtered rows.
- `select_groups(..., kind_input_id="eggnog_protein")` preserves `GroupId`.
- unmapped IDs are reported correctly.
- `write_tidy()` writes canonical `mapping.parquet`.
- gzip SQLite snapshots are decompressed under `dir_tmp`.
- invalid `kind_input_id` raises targeted `ValueError`.

## Real-Data Validation

Validated inputs:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/eggnog/mapper/5.0.2/raw/eggnog.db.gz
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/eggnog/cog/COG2024/raw/cog-24.fun.tab
```

Validated output:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/eggnog/mapper/5.0.2/tidy/mapping.parquet
```

Observed publication result on 2026-06-08:

- parquet size: about 949 MB
- row count: 118,683,777
- schema:
  - `EggnogProteinId`
  - `EggnogOgId`
  - `EggnogLevel`
  - `CogCategory`
  - `CogClass`
  - `CogName`
  - `OgDescription`

Observed runtime characteristics:

- compressed input `eggnog.db.gz`: about 6.4 GB
- temporary decompressed SQLite: about 41 GB
- temporary TSV peak: about 12 GB+
- temporary workspace peak under `/tmp/bioextract-eggnog-run`: about 57 GB
- process RSS observed around 2.8 GB
- thread count observed around 29

The temporary workspace was automatically cleaned after successful completion.
