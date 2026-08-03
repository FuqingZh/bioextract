# GO and KEGG Tidy Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract` owns local resource snapshot access. GO exposes a multi-relation
ontology and publishes it as DuckDB. KEGG BRITE and mapping products are
independent one-relation DuckDB capability profiles.

The implemented MVP covers:

- GO OBO to ontology tidy tables
- GO OBO subset membership and subset definition tables
- GO term selection by ID, namespace, and subset membership
- KEGG BRITE JSON to pathway tidy tables
- in-memory use through `build_tidy().frames`
- persisted GO use through `write_duckdb(path)`
- persisted KEGG use through `write_duckdb(path)`, `from_duckdb(path)`, and a
  fresh caller-owned read-only `connect()` connection

The MVP intentionally does not cover:

- GO GAF or GPAD annotation
- protein or gene identifier selection for GO
- ORA, GSEA, or clusterProfiler replacement logic
- projecting arbitrary GO annotations to GO slim ancestors

## Public API

```python
from bioextract import GODatabase, KEGGDatabase

go = GODatabase.from_obo("go-basic.obo")
go_tidy = go.build_tidy()
df_terms = go_tidy.frames["term"]
df_subsets = go.list_subsets()
df_goslim_generic = go.select_terms(
    subset_id="goslim_generic",
)
go_result = go.write_duckdb("out/go.duckdb")

kegg_tidy = KEGGDatabase.from_brite_json("br08901.json").build_tidy()
df_pathway = kegg_tidy.frames["pathway"]
kegg_result = KEGGDatabase.from_brite_json("br08901.json").write_duckdb(
    "out/kegg.duckdb"
)
published_kegg = KEGGDatabase.from_duckdb(kegg_result.path)
with published_kegg.connect() as connection:
    pathway_count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
```

Legacy directory writers remain only for migration. New callers use the
single-file writers above.

## Data Flow

`GODatabase` and `KEGGDatabase` are path-first resource handles. They validate
file existence during construction, but they do not parse the raw resource
until a read, selection, or write operation requires it.

GO OBO parsing is stanza-streamed through `scan_obo_term_records()`. The tidy
builder consumes records once into column buffers, then materializes Polars
frames at the artifact boundary.

The `edge` frame preserves parsed OBO parent and relationship edges. Derived
`ancestor_all` and `depth` frames use hierarchical relation types only
(`is_a`, `part_of`) so non-hierarchical relationships such as `has_part` do not
create graph cycles in subset OBO snapshots.

KEGG BRITE JSON currently uses the standard library JSON parser, so the JSON
tree is loaded before row traversal. The traversal and frame construction are
linear and deterministic. True streaming JSON parsing would require an
additional dependency such as `ijson`; that is out of scope for this MVP.

## Output Contract

GO ontology DuckDB tables:

```text
term
term_relation
term_synonym
term_xref
term_alternate_id
subset_membership
subset_definition
term_ancestor
term_depth
```

KEGG BRITE DuckDB tables:

```text
pathway
```

GO and KEGG provenance and table counts are stored in the metadata-v1
`_bioextract` relations. A reopened KEGG BRITE publication validates resource
identity, the `kegg-brite-json-v1` profile, the BRITE resource schema version,
the `brite` scope, and the exact `pathway` table/role/count before exposing
read-only SQL. Neither output requires a sidecar.
