# ReactomeDatabase Architecture

Version: v1.0
Date: 2026-07-30
Status: current

## Goal

`bioextract.reactome.ReactomeDatabase` provides path-first access to local Reactome
mapping snapshots. It supports annotation lookup and
standard enrichment inputs from Reactome open-data files without calling the
Reactome web services at runtime.

The implemented MVP covers:

- UniProt accession to Reactome pathway mapping
- Reactome pathway metadata
- Reactome pathway parent-child relations
- species-scoped extraction
- single-query and grouped protein selections
- `term2gene` and `term2name` frames for ORA and GSEA callers
- canonical single-DuckDB publication through the shared materialized-dataset
  contract
- composable input files, so annotation-only and metadata-only use cases do
  not need unrelated raw files

The MVP intentionally does not cover:

- Reactome online Analysis Service calls
- ReactomePA result compatibility
- enrichment p-value calculation
- reaction-level graph traversal
- cross-species orthology inference
- identifier mapping from gene symbol, Entrez, Ensembl, or STRING IDs

## Raw Inputs

The local v96 mapping snapshot currently has:

```text
UniProt2Reactome.txt
ReactomePathways.txt
ReactomePathwaysRelation.txt
```

`UniProt2Reactome.txt` is tab-separated with six columns:

```text
uniprot_id
reactome_pathway_id
reactome_url
pathway_name
evidence_code
species
```

`ReactomePathways.txt` is tab-separated with three columns:

```text
reactome_pathway_id
pathway_name
species
```

`ReactomePathwaysRelation.txt` is tab-separated with two columns:

```text
parent_reactome_pathway_id
child_reactome_pathway_id
```

The files are small enough for eager Polars reads in the first version. The
largest current input is `UniProt2Reactome.txt`, about 43 MB and 322,435 rows.

## Public API

```python
from bioextract import ReactomeDatabase

db = ReactomeDatabase.from_files(
    uniprot_mapping="UniProt2Reactome.txt",
    pathways="ReactomePathways.txt",
    relations="ReactomePathwaysRelation.txt",
)

selection = (
    db
    .with_species("Homo sapiens")
    .select_ids(["P04637", "Q9Y243", "MISSING"])
)

df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmatched_ids()

df_term2gene = db.with_species("Homo sapiens").extract_term2gene()
df_term2name = db.with_species("Homo sapiens").extract_term2name()
df_relations = db.with_species("Homo sapiens").extract_pathway_relations()
```

Grouped selections should mirror the existing STRINGdb and OmniPath shape:

```python
df_group_mapping = (
    db.with_species("Homo sapiens")
    .select_groups(
        {
            "TumorA": ["P04637", "Q9Y243"],
            "TumorB": ["P31749", "P42345"],
        }
    )
    .extract_mapping()
)
```

## Data Flow

`ReactomeDatabase.from_files()` accepts any non-empty combination of the three
raw files, validates provided file existence, then stores paths only. Parsing
happens when a frame is first needed.

Capability dependencies are explicit:

```text
extract_mapping              -> UniProt2Reactome.txt
extract_unmatched_ids   -> UniProt2Reactome.txt
extract_term2gene            -> UniProt2Reactome.txt
extract_term2name            -> ReactomePathways.txt
extract_pathway_relations    -> ReactomePathwaysRelation.txt
species-scoped relations     -> ReactomePathwaysRelation.txt + ReactomePathways.txt
```

If a caller invokes a capability whose backing file was not provided, the method
raises a targeted `ValueError`. This mirrors the STRINGdb behavior where aliases
and links can be supplied independently and missing-file failures occur at the
feature boundary.

`with_species(species)` returns a lightweight resource view with a species
filter. The filter should match the Reactome species display name exactly after
trimming whitespace. Case-insensitive convenience can be added later if there is
a demonstrated caller need, but exact matching keeps the first contract
auditable.

`select_ids(ids)` normalizes input IDs by trimming whitespace and dropping empty
values. The selected IDs are interpreted as UniProt accessions. The method does
not attempt gene-symbol or isoform conversion.

Extraction is DataFrame-first:

1. read the needed raw file once
2. apply species filter when present
3. join or filter with the normalized selected IDs
4. return deterministic, stable columns

The implementation may cache raw frames and extraction frames on the selection
object, following the existing STRINGdb and OmniPath pattern.

## Output Contract

`extract_mapping()` returns:

```text
input_id
uniprot_id
reactome_pathway_id
pathway_name
evidence_code
species
reactome_url
```

For grouped selections it prepends:

```text
group_id
input_id
uniprot_id
reactome_pathway_id
pathway_name
evidence_code
species
reactome_url
```

`extract_unmatched_ids()` returns:

```text
input_id
```

For grouped selections it returns:

```text
group_id
input_id
```

`extract_term2gene()` returns:

```text
reactome_pathway_id
uniprot_id
```

`extract_term2name()` returns:

```text
reactome_pathway_id
pathway_name
species
```

`extract_pathway_relations()` returns:

```text
parent_reactome_pathway_id
child_reactome_pathway_id
```

When species-scoped, relations should be limited to edges where both parent and
child exist in the species-scoped pathway metadata.

## Materialized Dataset

`build_tidy()` is optional but useful for resource publication and snapshot
inspection. Its internal frames are:

```text
mapping
pathway
relation
term2gene
term2name
```

With a partial snapshot, it builds only frames derivable from the provided
files. `term2gene` and `term2name` remain convenient in-memory enrichment
projections.

Suggested schema version:

```text
reactome-mapping-v0.1
```

Canonical publication uses `write_duckdb(path)`. `protein_pathway`,
`pathway`, and `pathway_relation` stay together; duplicate enrichment
projections are not stored. Provenance and row counts live in `_bioextract`.

`ReactomeDatabase.from_duckdb(path)` reopens full or partial publications that
conform to metadata v1 and `reactome-mapping-files-v1`. Reopen validation is
bounded to provenance, source capabilities, the exact table/view inventory,
and physical column schemas; recorded biological row counts are trusted rather
than recomputed. Existing selection, grouped-selection, enrichment, relation,
and species-scoping behavior is shared with source-backed handles.

Each `connect()` call on a reopened handle returns an independent native DuckDB
connection opened read-only and owned by the caller. The handle pins the
validated file identity and rejects access after atomic replacement. Source
handles do not expose native connections, and reopened handles do not republish
the validated database.

The pre-convergence verified v96 artifact is retained as a legacy baseline:

```text
reactome/mapping/v96/tidy/reactome.duckdb
```

It contains 322,435 `protein_pathway` rows, 23,498 `pathway` rows, and
23,612 `pathway_relation` rows. The prior multi-Parquet `tidy/reactome/`
directory is also a preserved migration artifact. Neither legacy path
overrides the current versioned CephFS convention, `tidy/data.duckdb`; formal
rebuild acceptance is pending.

## Why Not reactome2py

`reactome2py` is useful for online Reactome API calls. It is not a replacement
for this layer because `bioextract` needs deterministic local resource access,
offline operation, and snapshot-specific outputs. A future caller may offer a
separate online integration, but `ReactomeDatabase` should stay local-file-first.

## Implementation Notes

- Keep `ReactomeDatabase` under `src/bioextract/reactome/`.
- Export only `ReactomeDatabase`; abbreviated aliases are not supported.
- Keep parsing helpers module-level and prefix them with `read_`, `filter_`, or
  `extract_` according to their behavior.
- Do not introduce pandas or network dependencies.
- Prefer Polars expressions for filtering, joins, and deduplication.
