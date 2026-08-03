# EggNOGDatabase Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.eggnog.EggNOGDatabase` provides path-first access to local eggNOG mapper
resources and exposes stable selected protein-to-COG annotations for downstream
selection and enrichment-input extraction.

The current version covers:

- `eggnog.db.gz` or `eggnog.db` as the primary mapper database
- optional `cog-24.fun.tab` lookup for COG class and display name
- single and grouped selection by eggNOG protein ID

It intentionally does not cover:

- direct annotation from FASTA or DIAMOND hits
- emapper runtime execution
- KEGG/module/domain/pathway annotation beyond the local COG mapping path
- enrichment p-value calculation

## Raw Inputs

`EggNOGDatabase.from_sqlite(source, *, cog_functions=None, temp_dir=None)`
accepts an eggNOG mapper SQLite database or its gzip wrapper as `source`:

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
from bioextract import EggNOGDatabase

db = EggNOGDatabase.from_sqlite(
    "eggnog.db",
    cog_functions="cog-24.fun.tab",
)

selection = db.select_ids(
    ["9606.ENSP00000369497"],
)
df_selected = selection.extract_mapping()
df_unmapped = selection.extract_unmatched_ids()
```

Grouped selections prepend `GroupId` in the same style as other resource DBs.
The namespace is fixed to `eggnog_protein`.

Plain SQLite is queried directly and emits no warning. A gzip-wrapped source
emits one `UserWarning` at construction because it must be fully decompressed
for each uncached selection lookup. Decompression uses only temporary storage
under `temp_dir` when supplied, and that temporary content is removed after
successful and failed queries. Callers with repeated workloads should
decompress once themselves and pass the resulting `.db` file.

eggNOG remains direct-only: there is no database-level full mapping extractor,
publication writer, persistent cache, or public unpack helper. A selection
caches its mapping and unmatched outputs, so repeated extractors on that
selection do not repeat the SQLite lookup.

## Selected Output Contract

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

This keeps selected output auditable and easy to project into downstream term
tables without eagerly loading the complete snapshot.

## Selection Contract

Accepted `namespace` values:

```text
eggnog_protein
```

Single selection output prepends:

```text
InputId
InputNamespace
```

Grouped selection output prepends:

```text
GroupId
InputId
InputNamespace
```

`extract_unmatched_ids()` returns `InputId` for single selection and
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
[EggNOGDatabase test standard](../testing/eggnog-db.md) and the
[eggNOG 5.0.2 benchmark](../benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md).
