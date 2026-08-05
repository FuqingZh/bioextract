# ChEBIDatabase Architecture

Version: v1
Date: 2026-07-29
Status: current

## Purpose And Sources

`ChEBIDatabase` publishes a stable compound identity and ChEBI-owned
relationships from local official snapshots. FULL OBO is canonical. SDF only
adds molfile records, TSV remains a partial-build path, and the PostgreSQL dump
is an audit or future advanced-source input rather than a runtime dependency.
Optional ChemOnt OBO remains an independent `chemont_*` graph.

`from_obo(source, sdf=..., chemont_obo=...)` selects the ontology
representation and discovers exactly one OBO candidate from a directory or
archive; an SDF supplement may be discovered or explicitly replaced.
`from_table_files(source=None, ...)` selects the table representation and
merges explicit table-role replacements over a discovered directory/archive
profile. Compression and archive containers are detected internally, and a
mixed source never silently chooses between OBO and tables.

## Canonical Relations

The ChEBI publication contains `compound`, `secondary_id`, `compound_name`,
`compound_cross_reference`, `compound_relation`, `compound_structure`, and
`compound_wurcs`. It does not duplicate them as a `chebi_term*` family.
ChemOnt uses `chemont_term`, `chemont_term_relation`,
`chemont_term_synonym`, and `chemont_term_xref`.

All public ChEBI keys are complete `CHEBI:<number>` CURIE strings. Scalar
chemical representations such as SMILES, InChI, and InChIKey live on
`compound`; the potentially large molfile stays in `compound_structure`.
One-to-many names, xrefs, relations, structures, and WURCS values are extracted
separately so `extract_compounds()` remains one row per canonical compound.

## Domain And Native Access

`from_duckdb(path)` validates resource identity, metadata version, physical
CURIE typing, table inventory, and row counts. `select_compounds()` and
`select_groups()` resolve primary/secondary ChEBI IDs, exact InChI/InChIKey,
and dynamic external prefixes such as `kegg.compound` or `hmdb`.

Selections defer DuckDB work until an eager `extract_*()` terminal. They retain
caller and group lineage, support star/obsolete policy, return direct
relations, and traverse cycle-safe `is_a` ancestry or descendants.

`connect()` is the advanced escape hatch: every call returns a new native
`DuckDBPyConnection` opened with `read_only=True`. The library does not expose
thin `sql()`/`query()` wrappers, a shared connection, or a write option.

## Integrity And Publication

Missing or duplicate canonical ChEBI IDs fail the build with
`bioextract.errors.IntegrityError`. Orphan dependent records are skipped and
recorded as `foreign_key_violation` warnings in
`_bioextract.validation_issue`.

Metadata schema v1 has `metadata`, `source_file`, `table_info`,
`column_mapping`, and `validation_issue`. It records validation status and
issue count. Readers accept only this complete v1 shape, including the explicit
resource schema and source profile keys.

`write_duckdb(path, if_exists="fail")` builds a staging database, verifies it,
closes DuckDB, and atomically replaces the destination. The artifact has no
sidecar manifest.

## Non-goals

The adapter does not download releases, restore PostgreSQL, perform fuzzy name
or chemical similarity search, infer chemical relations, or merge ChEBI and
Rhea into an application-specific model.
