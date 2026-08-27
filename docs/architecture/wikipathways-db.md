# WikiPathwaysDatabase Architecture

Version: v1.1
Date: 2026-08-27
Status: current

## Goal

`bioextract.wikipathways.WikiPathwaysDatabase` provides source-first access to
local WikiPathways GMT snapshots. The first version is a local enrichment-input
layer: it reads one or more GMT files and emits pathway metadata through
`pathways()`, `pathway_genes()`, and `pathway_names()`.

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

Selected gene input is exact NCBI Entrez Gene text after surrounding whitespace
is trimmed. Empty values are dropped and exact duplicates collapse.
Pipe-bearing text is not treated as a UniProt representation and remains an
ordinary unmatched WikiPathways input.

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
from bioextract import WikiPathwaysDatabase

db = WikiPathwaysDatabase.from_gmt("wikipathways-20260510-gmt-*.gmt")
human = db.with_species("Homo sapiens")

lf_pathways = human.pathways()
lf_pathway_genes = human.pathway_genes()
lf_pathway_names = human.pathway_names()

selection = human.select_ids(["2687", "2678", "MISSING"])
lf_mapping = selection.mappings()
lf_unmapped = selection.unmatched_ids()
```

Canonical publications can be reopened without the original GMT files:

```python
published = WikiPathwaysDatabase.from_duckdb("tidy/wikipathways.duckdb")
with published.connect() as connection:
    pathway_count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
```

`from_duckdb()` requires metadata schema v2, the `wikipathways-gmt-v1` source
profile, and the current resource schema. It also requires official
content-derived release identity, the exact `pathway` and `pathway_gene`
canonical inventory, and their recorded roles, column provenance, and physical
schemas. Validation is bounded to metadata and catalogs: recorded biological
row counts must be non-negative but are not recounted during reopen.
`connect()` returns a fresh caller-owned native DuckDB connection in read-only
mode.

The reopened handle is pinned to the file identity that passed validation.
Every domain terminal, including a previously cached selection terminal,
rejects a vanished or atomically replaced publication and requires reopening.
Reopened tidy datasets cannot be republished because they do not own the GMT
source provenance needed to create a new canonical publication.

Grouped selections mirror the STRINGdb and Reactome style:

```python
df_mapping = (
    human.select_groups({"A": ["2687"], "B": ["435", "MISSING"]})
    .mappings()
)
df_mapping = df_mapping.collect()
```

## Source Resolution And Dataset Identity

`from_gmt(source, *, glob=True)` accepts one local string or
path-like path, a sequence of them, or glob expressions. `**` expressions
recurse into nested directories. With `glob=False`, each scalar or sequence
entry is literal. Empty sources, unmatched expressions, missing paths,
directories, and other non-files are rejected.

`with_species(species)` creates a view over the same source or reopened
publication. It applies the exact normalized species to `pathway`, then
semi-joins `pathway_gene` through the retained pathway IDs. This is the
species boundary that prevents an Entrez Gene ID shared by multiple species
from carrying a pathway across species.

Resolution is private to the WikiPathways package. Every match is normalized to
its actual physical file, duplicate physical files are rejected even when
reached through repeated inputs, overlapping expressions, or symlink aliases,
and the frozen file tuple is sorted deterministically before parsing. Filenames
and directories never supply release, species, Collection, Version, or schema
identity.

The complete unfiltered dataset must contain one official `Collection` and
globally unique `WikiPathwaysId` values. `Version` is derived from the
non-empty suffix of an official Collection field that starts with
`WikiPathways_`, and is also checked for one common value defensively.
Duplicate IDs within one file and across files are errors. One GMT may contain
several species.

The validated common Version is recorded as `release_version` with
`release_version_source="official_metadata"`. It comes only from parsed GMT
content, never from a filename, directory, or species-filtered result.

Construction resolves and freezes paths immediately. Content parsing and the
Collection, Version, and pathway-ID checks are deferred until the first
extraction, tidy build, or publication.

Every resolved GMT must independently contain at least one non-empty pathway
record. An empty or whitespace-only file is rejected with its path even when
other resolved files contain valid records.

## Output Contract

`pathways()` returns a native `polars.LazyFrame` with:

```text
wiki_pathways_id
pathway_name
species
collection
version
url
gene_count
```

`pathway_genes()` returns a native `polars.LazyFrame` with:

```text
wiki_pathways_id
gene_id
```

`pathway_names()` returns a native `polars.LazyFrame` with:

```text
wiki_pathways_id
pathway_name
species
collection
version
url
```

`mappings()` returns a native `polars.LazyFrame` with:

```text
input_id
gene_id
wiki_pathways_id
pathway_name
species
url
```

Grouped mapping prepends `group_id`.

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
Pathway names are not duplicated into another physical relation because the
canonical `pathway` relation already owns them.

## Species Filtering

`with_species()` is an exact row/content filter, not a source discovery or file
identity rule. When provided, pathway metadata is filtered by exact species
string, and `pathway_genes()` is filtered through the retained pathway IDs. Parsing,
Collection/Version validation, ID uniqueness checks, and provenance inventory
always cover every resolved file before this filter is applied.
