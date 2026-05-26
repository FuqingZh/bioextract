# GO and KEGG Tidy Architecture

## Goal

`bioextract` owns local resource snapshot access. The GO and KEGG tidy layer
turns raw local snapshots into stable Polars frames and optional parquet
artifacts. Manifest writing is optional.

The implemented MVP covers:

- GO OBO to ontology tidy tables
- KEGG BRITE JSON to pathway tidy tables
- in-memory use through `build_tidy().frames`
- persisted use through `build_tidy().write(dir_out)`
- optional manifest use through `build_tidy().write(dir_out, should_write_manifest=True)`

The MVP intentionally does not cover:

- GO GAF or GPAD annotation
- protein or gene identifier selection for GO
- ORA, GSEA, or clusterProfiler replacement logic
- GO slim mapping

## Public API

```python
from bioextract.go import GoDb
from bioextract.kegg import KeggDb

go_tidy = GoDb.from_obo("go-basic.obo").build_tidy()
df_terms = go_tidy.frames["term"]
go_report = go_tidy.write("out/go-basic")
go_report_with_manifest = go_tidy.write(
    "out/go-basic-archive",
    should_write_manifest=True,
)

kegg_tidy = KeggDb.from_brite_json("br08901.json").build_tidy()
df_pathway = kegg_tidy.frames["pathway"]
kegg_report = kegg_tidy.write("out/br08901")
```

`write_tidy(dir_out)` is available as a convenience wrapper around
`build_tidy().write(dir_out)`.

## Data Flow

`GoDb` and `KeggDb` are path-first resource handles. They validate file
existence and configured size limits during construction, but they do not parse
the raw resource until `build_tidy()`.

GO OBO parsing is stanza-streamed through `scan_obo_term_records()`. The tidy
builder consumes records once into column buffers, then materializes Polars
frames at the artifact boundary.

KEGG BRITE JSON currently uses the standard library JSON parser, so the JSON
tree is loaded before row traversal. The traversal and frame construction are
linear and deterministic. True streaming JSON parsing would require an
additional dependency such as `ijson`; that is out of scope for this MVP.

## Output Contract

GO ontology tidy output:

```text
term.parquet
edge.parquet
synonym.parquet
xref.parquet
alt_id.parquet
ancestor_all.parquet
depth.parquet
```

KEGG BRITE tidy output:

```text
pathway.parquet
```

When `should_write_manifest=True`, `manifest.json` is also written. The manifest
contains:

- `build_id`
- `schema_version`
- `generated_at`
- `sources`
- `assets`

Each asset entry contains:

- `path`
- `kind`
- `sha256`
- `row_count`
- `is_optional`

The output is intentionally flat because `bioextract` is library-first. The
prior `biotidy` `canonical/` and `derived/` directories are not part of the
new default contract.
