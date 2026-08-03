# KEGG Metabolic Database

Version: v0.1
Date: 2026-07-30
Status: current

## Purpose

`KEGGDatabase` publishes caller-supplied local KEGG compound, reaction,
enzyme, module, and global-link files as one provenance-aware DuckDB. This
metabolic product is independent of the KEGG BRITE and organism-mapping
Parquet products and does not change their schemas or release scopes.

## Construction And Publication

Use `from_metabolic_release()` for a complete release directory, its `raw`
directory, or a zip/tar archive. Use the explicit keyword-only
`from_metabolic_files()` roles for partial fixtures and nonstandard layouts.
Both constructors accept optional `release_version`. It records only a
caller-declared official release identity with
`release_version_source=caller`; directory, file, and archive names never
supply or validate release identity.
Entry batches are streamed record by record; excluded large fields are not
published. The writer stages locally beside the requested destination,
publishes metadata schema v1, verifies inventory and row counts read-only, and
then atomically commits `<resource>.duckdb`.

When a release is supplied as an archive, provenance records the original
archive display path and media type rather than an ephemeral extraction path.
Archive members remain subject to the shared safe-path and supported-format
checks.

The resource schema is `kegg-metabolic-v0.1`. Biological relations live in
`main`; `_bioextract` contains metadata, sources, table inventory, column
mapping, and validation issues. Partial inputs omit unavailable tables and
record their available source-role capabilities.

## Domain Semantics

Reaction equations retain their original text and publish ordered `left` and
`right` participants. KEGG reversibility is recorded without inventing
physiological substrate/product direction. Numeric coefficients have a numeric
projection; symbolic coefficients such as `n` and `(n+1)` remain text.

Module definitions publish a lossless ordered expression tree. Exact module
evaluation treats top-level sequence blocks as required, complexes as AND,
alternatives as OR, and optional nodes as satisfied. It returns counts,
completeness, and one-based missing block indexes without a caller-specific
threshold.

Selections accept the closed metabolic namespaces `kegg_compound`, `chebi`,
`pubchem`, `kegg_reaction`, `rhea`, `ec`, `ko`, `kegg_module`, and
`kegg_pathway`. They resolve canonical anchors and then select only directly
linked reactions. Extractors retain input/group lineage and project reaction
participants, compounds, ECs, KOs, modules, pathways, and cross-references.
Namespaces backed by cross-reference relations are available only when that
specific namespace has published values, not merely because the relation table
exists.

Real KEGG ENZYME obsolete entries are identified from
`ENTRY EC ... Obsolete Enzyme` together with `NAME` and `COMMENT`. Transferred
EC inputs recursively resolve replacement chains to active canonical ECs by
default; `include_obsolete=True` retains the exact historical EC. Deleted
entries and broken replacement targets remain distinguishable unmatched
outcomes, and missing replacement targets are persisted as validation issues.

## Native Access

`KEGGDatabase.from_duckdb(path)` validates identity, supported metadata and
resource schema versions, table inventory, and row-count parity.
`connect()` returns a new caller-owned native DuckDB connection opened with
`read_only=True`; arbitrary SQL remains a DuckDB concern.
