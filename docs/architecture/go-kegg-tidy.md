# GO and KEGG Tidy Architecture

Version: v1.2
Date: 2026-08-27
Status: current

## Goal

`bioextract` owns local resource snapshot access. GO exposes a multi-relation
ontology and publishes it as DuckDB. KEGG BRITE and mapping products are
independent one-relation DuckDB capability profiles.

The implemented MVP covers:

- GO OBO to ontology tidy tables
- GO OBO subset membership and subset definition tables
- GO term selection by ID, namespace, and subset membership
- GO ancestor selection with optional OBO subset projection and unmatched-ID
  accounting
- KEGG BRITE JSON to pathway tidy tables
- in-memory use through `build_tidy().frames`
- persisted GO use through `write_duckdb(path)`, `from_duckdb(path)`, and a
  fresh caller-owned read-only `connect()` connection
- persisted KEGG use through `write_duckdb(path)`, `from_duckdb(path)`, and a
  fresh caller-owned read-only `connect()` connection

The MVP intentionally does not cover:

- GO GAF or GPAD annotation
- protein or gene identifier selection for GO
- ORA, GSEA, or clusterProfiler replacement logic
- GO GAF/GPAD annotation projection and protein-to-GO mapping

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
ancestor_selection = go.select_ancestors(
    ["GO:0008150", "GO:1234567"],
    target_subset_id="goslim_generic",
    include_self=True,
)
lf_ancestors = ancestor_selection.ancestors()
lf_unmatched = ancestor_selection.unmatched_ids()
df_ancestors = lf_ancestors.collect()
df_unmatched = lf_unmatched.collect()
go_result = go.write_duckdb("out/go.duckdb")
published_go = GODatabase.from_duckdb(go_result.path)
with published_go.connect() as connection:
    term_count = connection.sql("SELECT count(*) FROM term").fetchone()[0]

kegg_tidy = KEGGDatabase.from_brite_json("br08901.json").build_tidy()
df_pathway = kegg_tidy.frames["pathway"]
kegg_result = KEGGDatabase.from_brite_json("br08901.json").write_duckdb(
    "out/kegg.duckdb"
)
published_kegg = KEGGDatabase.from_duckdb(kegg_result.path)
with published_kegg.connect() as connection:
    pathway_count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
```

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

Before graph derivation, a source build requires both endpoints of every
`is_a` and `part_of` edge to exist in the parsed term inventory. A dangling
edge raises contextual `ValueError` naming its child and parent. The equivalent
semantic defect in a reopened DuckDB is a publication-integrity failure and
`from_duckdb()` raises `IntegrityError`.

`select_terms()` and `list_subsets()` on a reopened publication execute their
resource-owned filters and aggregation in DuckDB. `select_ancestors()` also
pushes canonical/alternate ID resolution, hierarchical joins, optional self
rows, subset membership, and obsolete policy into DuckDB, then exposes the
bounded result through replayable Polars `LazyFrame` terminals. Source tidy
construction materializes the transitive ontology closure as the documented
global-context exception before wrapping its frames as replayable lazy
relations. `build_tidy()` remains the explicit complete-frame path. Source-
backed and publication-backed handles expose the same result schemas and
ordering even though their execution mechanisms differ.

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

The stable `build_tidy().frames` to physical-table mapping is intentional:

| Frame | DuckDB table |
| --- | --- |
| `term` | `term` |
| `edge` | `term_relation` |
| `synonym` | `term_synonym` |
| `xref` | `term_xref` |
| `alt_id` | `term_alternate_id` |
| `subset_membership` | `subset_membership` |
| `subset_definition` | `subset_definition` |
| `ancestor_all` | `term_ancestor` |
| `depth` | `term_depth` |

Both name sets are contract surfaces; changing either side requires a successor
resource schema rather than an incidental rename.

KEGG BRITE DuckDB tables:

```text
pathway
```

GO and KEGG provenance and table counts are stored in the metadata-v2
`_bioextract` relations. A reopened KEGG BRITE publication validates resource
identity, the `kegg-brite-json-v1` profile, the BRITE resource schema version,
the `brite` scope, and the exact `pathway` table/role/count before exposing
read-only SQL. Neither output requires a sidecar.

`from_duckdb()` captures the KEGG BRITE file identity before profile
inspection and confirms the same identity after full validation. The reopened
handle checks that identity around each later `connect()` open. Replacing or
removing the path invalidates the old handle and requires an explicit
`from_duckdb()` reopen; it never silently changes the publication represented
by the handle.

A reopened GO publication validates the GO resource identity,
`gene-ontology-obo-v1` profile, ontology schema version, exact nine-table
ontology capability, semantic roles, and physical column schemas. Validation
uses the bounded catalog and provenance relations and does not recount the
biological tables. The handle pins the validated file identity, so replacing
the path requires an explicit reopen before domain reads or native SQL access.
