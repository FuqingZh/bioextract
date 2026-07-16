# EggnogDb Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.eggnog.EggnogDb` provides path-first access to local eggNOG mapper
resources and exposes a stable protein-to-COG annotation table for downstream
selection, enrichment-input extraction, and resource publication.

The current version covers:

- `eggnog.db.gz` or `eggnog.db` as the primary mapper database
- optional `cog-24.fun.tab` lookup for COG class and display name
- full mapping extraction
- single and grouped selection by eggNOG protein ID
- flat tidy writing to `mapping.parquet`

It intentionally does not cover:

- direct annotation from FASTA or DIAMOND hits
- emapper runtime execution
- KEGG/module/domain/pathway annotation beyond the local COG mapping path
- enrichment p-value calculation

## Raw Inputs

`EggnogDb.from_files()` accepts exact file paths:

```text
eggnog.db.gz
cog-24.fun.tab
```

`cog-24.fun.tab` is optional. Without it, the output still keeps:

- `EggnogProteinId`
- `EggnogOgId`
- `EggnogLevel`
- `CogCategory`
- `OgDescription`

and leaves `CogClass` / `CogName` null.

## Public API

```python
from bioextract.eggnog import EggnogDb

db = EggnogDb.from_files(
    file_eggnog_db="eggnog.db.gz",
    file_cog_fun="cog-24.fun.tab",
)

df_mapping = db.extract_mapping()

selection = db.select_ids(
    ["9606.ENSP00000369497"],
    kind_input_id="eggnog_protein",
)
df_selected = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()
```

Grouped selections prepend `GroupId` in the same style as other resource DBs.

## Output Contract

The wide mapping table exposes:

```text
EggnogProteinId
EggnogOgId
EggnogLevel
CogCategory
CogClass
CogName
OgDescription
```

Many-to-many expansion is preserved:

- one protein may map to multiple OGs
- one OG may emit multiple `CogCategory` rows

This keeps the output auditable and easy to project into downstream term tables.

## Tidy Writing

`write_tidy()` emits:

```text
mapping.parquet
manifest.json
```

The writer does not materialize the full table in Python memory before parquet
write. The implemented path is:

1. open local eggNOG SQLite, decompressing `.gz` to `dir_tmp` when needed
2. stream expanded rows into a temporary TSV
3. scan the TSV lazily with Polars
4. `sink_parquet()` into the final tidy artifact

`build_tidy()` is intentionally not offered as a stable full-resource lazy
dataset handle for the SQLite source, because its lazy plan depends on a
materialized intermediate TSV with temporary lifetime.

## Selection Contract

Accepted `kind_input_id` values:

```text
eggnog_protein
```

Single selection output prepends:

```text
InputId
KindInputId
```

Grouped selection output prepends:

```text
GroupId
InputId
KindInputId
```

`extract_unmapped_input_ids()` returns `InputId` for single selection and
`GroupId, InputId` for grouped selection.

## Real Snapshot Status

The local mapper snapshot currently used in validation is:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/eggnog/mapper/5.0.2/raw/eggnog.db.gz
```

The local COG function lookup is:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/eggnog/cog/COG2024/raw/cog-24.fun.tab
```

Validation details are recorded in the
[EggnogDb test standard](../testing/eggnog-db.md) and the
[eggNOG 5.0.2 benchmark](../benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md).
