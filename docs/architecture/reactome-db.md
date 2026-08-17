# ReactomeDatabase Architecture

Version: v1.2
Date: 2026-08-18
Status: current

## Goal

`bioextract.reactome.ReactomeDatabase` provides offline, path-first access to
caller-supplied Reactome snapshots. It parses the official local files,
exposes replayable Polars relations, and can publish or reopen one
provenance-aware metadata-v2 DuckDB. It does not download data, call Reactome
web services, or calculate enrichment statistics.

The current implemented capability set is `reactome-mapping-v0.5`:

- twelve six-column identifier-to-Reactome mapping roles across `uniprot`,
  `ncbi`, `chebi`, and `gtop`;
- `pathway` and `pathway_relation` metadata roles;
- human `complex_pathway` and `ewas_pathway` entity membership; and
- human `pathway_gene_set` membership from the official one-member GMT zip.

Every source role has its own physical table. There is no union table, target
column, compatibility alias, or filename-derived schema identity.

## Construction And Roles

`from_files()` accepts explicit keyword-only roles. The mapping arguments are:

```text
uniprot_mapping       UniProt2Reactome.txt
uniprot_all_levels    UniProt2Reactome_All_Levels.txt
uniprot_reactions     UniProt2ReactomeReactions.txt
ncbi_mapping          NCBI2Reactome.txt
ncbi_all_levels       NCBI2Reactome_All_Levels.txt
ncbi_reactions        NCBI2ReactomeReactions.txt
chebi_mapping         ChEBI2Reactome.txt
chebi_all_levels      ChEBI2Reactome_All_Levels.txt
chebi_reactions       ChEBI2ReactomeReactions.txt
gtop_mapping          GtoP2Reactome.txt
gtop_all_levels       GtoP2Reactome_All_Levels.txt
gtop_reactions        GtoP2ReactomeReactions.txt
```

The remaining roles are `pathways`, `relations`, `complex_pathways`,
`ewas_pathways`, and `pathway_gene_sets`. Any non-empty combination is valid;
missing roles remain unavailable capabilities. Every supplied path must be a
regular file. Release identity is optional caller input and is never inferred
from a path, basename, directory, or manifest.

## Mapping Contract

All twelve mapping sources are literal, headerless six-field TSV files. CSV
quoting is disabled, ragged or empty records fail closed, exact six-field
duplicates are removed, and rows differing in evidence remain distinct. All
canonical values are strings.

Pathway mapping tables use this source-specific schema and order:

```text
<namespace>_id
reactome_pathway_id
reactome_url
pathway_name
evidence_code
species
```

Reaction mapping tables use:

```text
<namespace>_id
reactome_reaction_id
reactome_url
reaction_name
evidence_code
species
```

The public whole-resource relations are:

```python
db.pathway_mappings(namespace="uniprot", pathway_level="lowest_level")
db.reaction_mappings(namespace="uniprot")
db.pathway_genes(pathway_level="lowest_level")  # UniProt only
```

`pathway_mappings()` never unions namespaces or levels. `reaction_mappings()`
is evidence of identifier-to-event membership only; it does not claim reaction
participants, direction, catalysts, regulation, or topology.

Selection dimensions are explicit:

```python
db.select_ids(ids, namespace="uniprot", target="pathway", pathway_level=None)
db.select_groups(
    ids_by_group, namespace="uniprot", target="pathway", pathway_level=None
)
```

For pathway targets, `None` resolves to `lowest_level`; `all_levels` remains
explicit. Reaction targets require `pathway_level=None`. Selected output keeps
the source-specific identifier column. Existing UniProt pathway selection
retains its established column order; reaction selections use
`input_id, <namespace>_id, reactome_reaction_id, reaction_name, evidence_code,
species, reactome_url`.

Input normalization is namespace-aware:

| Namespace | Accepted input | Lookup | Public `input_id` |
| --- | --- | --- | --- |
| `uniprot` | accession or pipe form | accession | accession |
| `ncbi` | any non-empty trimmed official text | same text | same text |
| `chebi` | decimal or `CHEBI:<digits>` | decimal text | `CHEBI:<integer>` |
| `gtop` | decimal text | same text | same text |

Invalid non-empty ChEBI/GtoP inputs raise `ValueError`; they are not reported
as unmatched. NCBI identifiers are deliberately not restricted to decimal Gene
IDs, because the official snapshot also contains GenBank/RefSeq-style values.

## Human Relations And GMT

`complex_pathways()` returns
`reactome_complex_id, reactome_pathway_id, top_level_reactome_pathway_id`.
`ewas_pathways()` returns the analogous `reactome_ewas_id` relation. The
headered source columns are exact and are recorded in
`_bioextract.column_mapping`; no entity-prefix species inference is used.

`pathway_gene_sets()` returns
`reactome_pathway_id, gene_set_name, gene_symbol` from exactly one regular
`ReactomePathways.gmt` member in the supplied zip. The archive is inspected
without extraction; encrypted, nested, parent-traversing, extra-member,
invalid-UTF-8, CRC, and decompression failures fail closed. GMT labels and gene
symbols are preserved as opaque source tokens. A pathway cannot have multiple
distinct GMT labels in one source.

Complex, EWAS, and GMT are file-scoped human relations: unscoped and
`with_species("Homo sapiens")` return all rows, while any other species returns
an empty relation with the stable schema. They have no selected/grouped API or
per-row species column.

## Validation And Publication

When pathway metadata is present, mapping rows, entity endpoints, and GMT
pathway IDs missing from human metadata are preserved and produce visible
`missing_pathway_metadata` warnings. Mapping closure is checked per namespace
when lowest-level, all-level, and hierarchy roles are all available. Entity
top-level claims are checked for reflexive ancestry when hierarchy is present;
cycles or contradictory claims are fatal. Absent comparison roles make the
check unavailable rather than synthesizing data.

The v0.5 publication identity is:

```text
bioextract.metadata_schema_version = 2
bioextract.source_schema_profile   = reactome-mapping-files-v5
bioextract.resource_schema_version = reactome-mapping-v0.5
```

`write_duckdb()` publishes exactly the available biological roles plus the five
shared `_bioextract` metadata tables. `_bioextract.source_file` is the
capability inventory. The GMT source is recorded as `application/zip`; all
other Reactome sources are `text/tab-separated-values`. Reopen validates exact
role/table presence, ordered physical schemas, expected entity column lineage,
media types, validation-state consistency, and the pinned file identity without
recounting biological rows. v0.4 and earlier Reactome publications are
intentionally rejected.

## Formal Boundary

Source parsing and temporary publications are separate from formal resource
replacement. A formal v96 delivery must build outside the live resource path,
validate the complete role set, preserve a rollback artifact, replace the
target atomically, and regenerate the release catalog from the exact package
and publication. Package release, catalog admission, downstream activation,
and deletion of old artifacts remain separately controlled operations.

## Implementation Notes

- Keep `ReactomeDatabase` under `src/bioextract/reactome/` and export only the
  database class.
- Keep the role registry and parsers private to Reactome; do not add a generic
  query facade or network dependency.
- Use Polars relations for source/publication query execution and deterministic
  canonical ordering.
