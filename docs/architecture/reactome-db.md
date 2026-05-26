# ReactomeDb Architecture

## Goal

`bioextract.reactome.ReactomeDb` provides path-first access to local Reactome
mapping snapshots. The first version should support annotation lookup and
standard enrichment inputs from Reactome open-data files without calling the
Reactome web services at runtime.

The implemented MVP should cover:

- UniProt accession to Reactome pathway mapping
- Reactome pathway metadata
- Reactome pathway parent-child relations
- species-scoped extraction
- single-query and grouped protein selections
- `term2gene` and `term2name` frames for ORA and GSEA callers
- optional tidy parquet writing through the shared `TidyDataset` contract
- composable input files, so annotation-only and metadata-only use cases do
  not need unrelated raw files

The MVP intentionally should not cover:

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
from bioextract.reactome import ReactomeDb

db = ReactomeDb.from_files(
    file_uniprot2reactome="UniProt2Reactome.txt",
    file_pathways="ReactomePathways.txt",
    file_relations="ReactomePathwaysRelation.txt",
)

selection = (
    db
    .with_species("Homo sapiens")
    .select_ids(["P04637", "Q9Y243", "MISSING"])
)

df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()

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

`ReactomeDb.from_files()` accepts any non-empty combination of the three raw
files, validates provided file existence and configured size limits, then stores
paths only. Parsing happens when a frame is first needed.

Capability dependencies are explicit:

```text
extract_mapping              -> UniProt2Reactome.txt
extract_unmapped_input_ids   -> UniProt2Reactome.txt
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
InputId
UniProtId
ReactomePathwayId
PathwayName
EvidenceCode
Species
ReactomeUrl
```

For grouped selections it prepends:

```text
GroupId
InputId
UniProtId
ReactomePathwayId
PathwayName
EvidenceCode
Species
ReactomeUrl
```

`extract_unmapped_input_ids()` returns:

```text
InputId
```

For grouped selections it returns:

```text
GroupId
InputId
```

`extract_term2gene()` returns:

```text
ReactomePathwayId
UniProtId
```

`extract_term2name()` returns:

```text
ReactomePathwayId
PathwayName
Species
```

`extract_pathway_relations()` returns:

```text
ParentReactomePathwayId
ChildReactomePathwayId
```

When species-scoped, relations should be limited to edges where both parent and
child exist in the species-scoped pathway metadata.

## Tidy Dataset

`build_tidy()` should be optional but useful for resource publication and
snapshot inspection. With all files present it should produce:

```text
mapping.parquet
pathway.parquet
relation.parquet
term2gene.parquet
term2name.parquet
```

With a partial snapshot, it should emit only frames that can be derived from the
provided files. For example, a mapping-only snapshot emits `mapping.parquet` and
`term2gene.parquet`; a pathways-only snapshot emits `pathway.parquet` and
`term2name.parquet`.

Suggested schema version:

```text
reactome-mapping-v0.1
```

Like GO and KEGG, the tidy writer should use flat output paths and optional
`manifest.json` through `TidyDataset.write()`.

## Why Not reactome2py

`reactome2py` is useful for online Reactome API calls. It is not a replacement
for this layer because `bioextract` needs deterministic local resource access,
offline operation, and snapshot-specific outputs. A future caller may offer a
separate online integration, but `ReactomeDb` should stay local-file-first.

## Implementation Notes

- Keep `ReactomeDb` under `src/bioextract/reactome/`.
- Export only `ReactomeDb`, `ReactomeResourceLimits`, and `ReactomeTidyDataset`
  from the package root.
- Keep parsing helpers module-level and prefix them with `read_`, `filter_`, or
  `extract_` according to their behavior.
- Do not introduce pandas or network dependencies.
- Prefer Polars expressions for filtering, joins, and deduplication.
