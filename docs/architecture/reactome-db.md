# ReactomeDatabase Architecture

Version: v1.1
Date: 2026-08-17
Status: current

## Goal

`bioextract.reactome.ReactomeDatabase` provides path-first access to caller-
supplied local Reactome mapping snapshots. It parses official TSV files,
exposes native Polars relations for annotation and enrichment callers, and can
publish or reopen one provenance-aware DuckDB. It does not call Reactome web
services or calculate enrichment statistics.

P1 implements two explicit UniProt-to-pathway capabilities:

- `lowest_level`: `UniProt2Reactome.txt`, the preserved default;
- `all_levels`: `UniProt2Reactome_All_Levels.txt`, the explicit hierarchy-
  expanded relation.

The two capabilities never share an ambiguous union table. Reaction,
NCBI Gene, ChEBI, GtoP, Complex, EWAS, and GMT roles remain deferred.

## Raw Inputs

The constructor accepts any non-empty combination of these caller-declared
roles:

```text
uniprot_mapping     -> UniProt2Reactome.txt
uniprot_all_levels  -> UniProt2Reactome_All_Levels.txt
pathways            -> ReactomePathways.txt
relations           -> ReactomePathwaysRelation.txt
```

Both UniProt mapping files are literal six-column TSV relations. Their columns
are preserved in this order:

```text
uniprot_id
reactome_pathway_id
reactome_url
pathway_name
evidence_code
species
```

The shared mapping-family parser disables CSV quoting, rejects ragged records
and empty required fields, retains literal quote characters, removes only
exact duplicate six-column rows, and retains rows that differ by
`evidence_code`.

`release_version` is optional caller-declared identity. It is trimmed and must
be non-empty when supplied. It is never inferred from a directory, filename,
manifest, or modification time.

## Public API

```python
from bioextract import ReactomeDatabase

db = ReactomeDatabase.from_files(
    uniprot_mapping="UniProt2Reactome.txt",
    uniprot_all_levels="UniProt2Reactome_All_Levels.txt",
    pathways="ReactomePathways.txt",
    relations="ReactomePathwaysRelation.txt",
    release_version="96",
)

lowest = db.pathway_mappings()
all_levels = db.pathway_mappings(pathway_level="all_levels")
lowest_genes = db.pathway_genes()
all_level_genes = db.pathway_genes(pathway_level="all_levels")

selection = db.with_species("Homo sapiens").select_ids(
    ["P04637", "MISSING"],
    pathway_level="all_levels",
)
lf_mapping = selection.mappings()
lf_unmapped = selection.unmatched_ids()
lf_names = db.with_species("Homo sapiens").pathway_names()
lf_relations = db.with_species("Homo sapiens").pathway_relations()
```

`namespace` and `target` are explicit selection dimensions on
`pathway_mappings()`, `select_ids()`, and `select_groups()`. P1 implements
only `namespace="uniprot"` and `target="pathway"`; invalid values raise
`ValueError`. A valid level whose source or publication table is absent raises
`ValueError` for a source handle or `CapabilityError` for a publication
handle.

Calls that omit `pathway_level` remain lowest-level UniProt pathway calls.
There is no implicit hierarchy closure and no silent union across levels.

`with_species()` applies an exact trimmed Reactome species display-name filter
before selection or enrichment deduplication. Pathway relations in a species
scope retain only edges whose two endpoints exist in that species' pathway
metadata.

## Output Contract

`pathway_mappings()` returns the six canonical mapping columns in source order:

```text
uniprot_id
reactome_pathway_id
reactome_url
pathway_name
evidence_code
species
```

`pathway_genes()` returns distinct, sorted pairs:

```text
reactome_pathway_id
uniprot_id
```

Canonical mapping relations retain evidence-distinct rows. Only the explicit
`pathway_genes()` projection deduplicates to a pathway/UniProt pair.

Selection mappings retain the existing shape:

```text
input_id
uniprot_id
reactome_pathway_id
pathway_name
evidence_code
species
reactome_url
```

Grouped mappings prepend `group_id`; unmatched outputs remain `input_id` or
`group_id, input_id`.

## Data Flow And Capability Boundaries

`from_files()` validates that every supplied path is a regular file and stores
paths only. Parsing is performed by the capability that needs the relation.
The all-level path is never used to satisfy the lowest-level role, and a
missing all-level source is not synthesized from `pathway_relation`.

The role-to-table boundary is private to the Reactome implementation:

```text
uniprot_pathway_lowest_level -> UniProt2Reactome.txt
uniprot_pathway_all_level    -> UniProt2Reactome_All_Levels.txt
pathway                      -> ReactomePathways.txt
pathway_relation             -> ReactomePathwaysRelation.txt
```

Source-backed and reopened handles use the same public query semantics. A
reopened handle validates metadata, source-role inventory, exact physical
tables, and column schemas without recounting biological tables. Each
`connect()` call returns an independent native DuckDB connection opened
read-only and pins the validated file identity.

## Materialized Dataset

`build_tidy()` exposes only the canonical role frames that are available:

```text
uniprot_pathway_lowest_level
uniprot_pathway_all_level
pathway
pathway_relation
```

`pathway_genes()` and `pathway_names()` are public lazy projections, not
redundant publication relations. `write_duckdb()` publishes exactly the available
canonical roles and the shared five `_bioextract` metadata-v2 relations.

The P1 publication identity is:

```text
bioextract.metadata_schema_version = 2
bioextract.source_schema_profile   = reactome-mapping-files-v2
bioextract.resource_schema_version = reactome-mapping-v0.2
```

Mapping rows whose pathway IDs are absent from supplied pathway metadata remain
published. When the affected mapping and `pathway` roles are both present,
publication records one `missing_pathway_metadata` warning per distinct
`(mapping_role, reactome_pathway_id)` in
`_bioextract.validation_issue`.

When lowest-level, all-level, and pathway-relation roles are all present, the
publisher compares all-level keys with the reflexive hierarchy closure of the
lowest-level keys at
`(uniprot_id, reactome_pathway_id, evidence_code, species)` grain. Any
mismatch is fatal; absence of a comparison role is simply an unavailable
check.

The old `protein_pathway` table, metadata-v1 reader, v0.1 reader, compatibility
view, and automatic migration are intentionally absent before v1.0.

## Formal Publication Boundary

P1 code and temporary local publications do not replace the formal v96
CephFS artifact. A later release must build outside the formal resource path,
validate all-level closure and warning preservation, then stage and replace an
artifact atomically under its own release authority. Package release, catalog
admission, downstream activation, and old-artifact deletion are separate
decisions.

## Why Not reactome2py

`reactome2py` is useful for online Reactome API calls. It is not a replacement
for this layer because `bioextract` needs deterministic local access,
offline operation, and snapshot-specific outputs.

## Implementation Notes

- Keep `ReactomeDatabase` under `src/bioextract/reactome/`.
- Export only `ReactomeDatabase`; selection classes remain implementation
  types.
- Keep mapping-family parsing shared and module-local; do not add a generic
  query facade or network dependency.
- Prefer Polars expressions for filtering, joining, closure validation, and
  deduplication.
