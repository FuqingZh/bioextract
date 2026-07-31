# WikiPathwaysDatabase Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.wikipathways.WikiPathwaysDatabase` provides source-first access to
local WikiPathways GMT snapshots. The first version is a local enrichment-input
layer: it reads one or more GMT files and emits pathway metadata plus
`term2gene`/`term2name` frames.

The MVP covers:

- single- and multiple-species WikiPathways GMT files
- NCBI Entrez Gene ID gene sets
- pathway metadata extraction
- single-query and grouped gene selections
- unmapped input ID reporting
- multi-relation DuckDB publication

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
from bioextract.wikipathways import WikiPathwaysDatabase

db = WikiPathwaysDatabase.from_gmt(
    "wikipathways-20260510-gmt-*.gmt",
    species="Homo sapiens",
)

df_pathway = db.extract_pathway()
df_term2gene = db.extract_term2gene()
df_term2name = db.extract_term2name()

selection = db.select_ids(["2687", "2678", "MISSING"])
df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmatched_ids()
```

Grouped selections mirror the STRINGdb and Reactome style:

```python
df_mapping = (
    db.select_groups({"A": ["2687"], "B": ["435", "MISSING"]})
    .extract_mapping()
)
```

## Source Resolution And Dataset Identity

`from_gmt(source, *, species=None, glob=True)` accepts one local string or
path-like path, a sequence of them, or glob expressions. `**` expressions
recurse into nested directories. With `glob=False`, each scalar or sequence
entry is literal. Empty sources, unmatched expressions, missing paths,
directories, and other non-files are rejected.

Resolution is private to the WikiPathways package. Every match is normalized to
its actual physical file, duplicate physical files are rejected even when
reached through repeated inputs, overlapping expressions, or symlink aliases,
and the frozen file tuple is sorted deterministically before parsing. Filenames
and directories never supply release, species, Collection, Version, or schema
identity.

The complete unfiltered dataset must contain one official `Collection` and
globally unique `WikiPathwaysId` values. `Version` is derived from the
`WikiPathways_` suffix of the official Collection field and is also checked for
one common value defensively. Duplicate IDs within one file and across files
are errors. One GMT may contain several species.

The validated common Version is recorded as `release_version` with
`release_version_source="official_metadata"`. It comes only from parsed GMT
content, never from a filename, directory, or species-filtered result.

Construction resolves and freezes paths immediately. Content parsing and the
Collection, Version, and pathway-ID checks are deferred until the first
extraction, tidy build, or publication.

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

## Publication

`write_duckdb(path)` publishes related pathway and membership relations:

```text
pathway
pathway_gene
```

Schema version:

```text
wikipathways-gmt-v0.1
```

Every resolved actual file is stored in `_bioextract` under a deterministic
unique logical source name, including files that contribute no rows after
filtering. The content-derived release Version, table roles, and row counts are
stored there as before.
`term2name` is not duplicated because the canonical `pathway` relation already
owns pathway names.

## Species Filtering

`species=` is an exact row/content filter, not a source discovery or file
identity rule. When provided, pathway metadata is filtered by exact species
string, and `term2gene` is filtered through the retained pathway IDs. Parsing,
Collection/Version validation, ID uniqueness checks, and provenance inventory
always cover every resolved file before this filter is applied.
