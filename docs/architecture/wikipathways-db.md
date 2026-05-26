# WikiPathwaysDb Architecture

## Goal

`bioextract.wikipathways.WikiPathwaysDb` provides path-first access to local
WikiPathways GMT snapshots. The first version is a local enrichment-input layer:
it reads species-specific GMT files and emits pathway metadata plus
`term2gene`/`term2name` frames.

The MVP covers:

- species-specific WikiPathways GMT files
- NCBI Entrez Gene ID gene sets
- pathway metadata extraction
- single-query and grouped gene selections
- unmapped input ID reporting
- flat tidy parquet writing

It intentionally does not cover:

- online WikiPathways API calls
- GPML graph parsing
- pathway topology or interaction extraction
- gene symbol, Ensembl, UniProt, or STRING ID conversion
- enrichment p-value calculation

## Official Format Contract

WikiPathways publishes monthly GMT files under:

```text
https://data.wikipathways.org/current/gmt/
```

The rWikiPathways `readPathwayGMT()` reference describes the returned data as
pathway-gene associations. The pathway fields include pathway name, version,
identifier, and organism, while gene content is provided as NCBI Entrez Gene
identifiers.

The local 2026-05-10 files follow this GMT layout:

```text
PathwayName%Collection%WikiPathwaysId%Species<TAB>Url<TAB>EntrezGeneId...
```

Example:

```text
Glutathione metabolism%WikiPathways_20260510%WP100%Homo sapiens	https://www.wikipathways.org/instance/WP100	2687	2678
```

## Public API

```python
from bioextract.wikipathways import WikiPathwaysDb

db = WikiPathwaysDb.from_gmt(
    "wikipathways-20260510-gmt-Homo_sapiens.gmt",
    species="Homo sapiens",
)

df_pathway = db.extract_pathway()
df_term2gene = db.extract_term2gene()
df_term2name = db.extract_term2name()

selection = db.select_ids(["2687", "2678", "MISSING"])
df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()
```

Grouped selections mirror the STRINGdb and Reactome style:

```python
df_mapping = (
    db.select_groups({"A": ["2687"], "B": ["435", "MISSING"]})
    .extract_mapping()
)
```

## Output Contract

`extract_pathway()` returns:

```text
WikiPathwaysId
PathwayName
Species
Collection
Version
Url
GeneCount
```

`extract_term2gene()` returns:

```text
WikiPathwaysId
GeneId
```

`extract_term2name()` returns:

```text
WikiPathwaysId
PathwayName
Species
Collection
Version
Url
```

`extract_mapping()` returns:

```text
InputId
GeneId
WikiPathwaysId
PathwayName
Species
Url
```

Grouped mapping prepends `GroupId`.

## Tidy Dataset

`build_tidy()` emits flat parquet assets:

```text
pathway.parquet
term2gene.parquet
term2name.parquet
```

Schema version:

```text
wikipathways-gmt-v0.1
```

No `canonical/` or `derived/` subdirectories are written. Asset `kind` is kept
only in `manifest.json` metadata.

## Species Filtering

GMT files are expected to be species-specific, but `species=` is still accepted
as a guard. When provided, pathway metadata is filtered by exact species string,
and `term2gene` is filtered through the retained pathway IDs.
