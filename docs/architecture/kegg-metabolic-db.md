# KEGG Metabolic Database

Version: v0.1
Date: 2026-08-07
Status: current

## Purpose

`KEGGDatabase` publishes caller-supplied local KEGG compound, reaction,
enzyme, module, and global-link files as one provenance-aware DuckDB. This
metabolic product is independent of the one-relation KEGG BRITE and
organism-mapping DuckDB products and does not change their schemas or release
scopes.

## Construction And Publication

Use `from_metabolic_files(source=None, ...)`. Without `source`, explicit roles
retain partial-profile behavior and missing roles become absent capabilities.
With `source`, the constructor accepts one release directory, its `raw`
directory, a parent containing one nested `raw` layout, or a zip/tar archive.
Zero plausible layouts reject as incomplete and multiple layouts reject as
ambiguous.

Every explicit non-`None` role replaces the complete matching discovery;
`None` alone permits discovery. Entry roles accept one batch, a directory, or
a sequence. Scalar list and relation roles accept one file or a directory.
Every explicit value must resolve to at least one file, and the final inventory
must not reuse one physical file across roles. Source-backed construction
requires all four entry collections and every relation role after overlays.
The four `*_list` roles are optional; list/entry parity is checked when a list
is present.

The constructor accepts optional `release_version`. It records only a
caller-declared official release identity with
`release_version_source=caller`; directory, file, and archive names never
supply or validate release identity.
Entry batches are streamed record by record; excluded large fields are not
published. The writer stages locally beside the requested destination,
publishes metadata schema v1, verifies inventory and row counts read-only, and
then atomically commits the requested DuckDB destination.

Directory-backed provenance contains only the final retained or replacement
files. When an archive contributes at least one retained role, provenance
records its original display path and each explicit overlay rather than
ephemeral extraction paths. A fully replaced archive is omitted. Archive
members remain subject to safe-path and supported-format checks.

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
