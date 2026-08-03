# Materialized Dataset Contract

Version: v1.0
Date: 2026-07-29
Status: current

## Purpose

This contract defines how `bioextract` publishes validated local snapshots. It
is subordinate to the repository-wide
[Domain Access Architecture](20260729-v1.0-domain-access-architecture.md):
materialization is an execution strategy, not the product purpose.

## Output Unit And Strategies

There are exactly two storage strategies:

- official/native direct access when the upstream representation is fit; or
- one bioextract-owned DuckDB per materialized logical product, regardless of
  relation count.

A partial capability of a multi-relation resource retains the same DuckDB
container. Missing inputs produce absent tables, not a different packaging
model and not misleading empty observations.

Resource publication writers use:

```python
file_out = "tidy/data.duckdb"
result = source.write_duckdb(file_out)
```

Here `file_out` is the positional destination value; it is not a promised
keyword name. Readers use `XDatabase.from_duckdb(path)`, and every `connect()`
call returns a fresh caller-owned connection opened read-only. The library does
not infer a semantic filename. The versioned CephFS convention is
`tidy/data.duckdb`, but filenames are not schema identity.

## Naming

Tables, views, and generated fields use singular lowercase `snake_case` and do
not contain hyphens.

Fields use a source-first rule:

- one-to-one columns copied from an official two-dimensional table retain the
  official header;
- fields parsed from OBO, RDF, XML, SDF, headerless inputs, joins, prefix
  decomposition, or other derivations use stable `snake_case`;
- source and normalized duplicate columns are not both published;
- empty, duplicate, or case-insensitively conflicting headers receive the
  smallest deterministic query-safe mapping;
- every necessary mapping is recorded in provenance.

An upstream official header change is therefore a schema change, not silently
hidden by a broad normalization layer.

## DuckDB provenance

Biological relations live in `main`. The `_bioextract` schema is an
application-owned internal provenance namespace, not a DuckDB or SQL reserved
name:

```text
main.<domain_table>

_bioextract.metadata
_bioextract.source_file
_bioextract.table_info
_bioextract.column_mapping
_bioextract.validation_issue
```

`metadata` stores resource identity, resource schema version, metadata schema
version, scope, package version, and generation time. `source_file` stores
logical source name, display path, bytes, media type, and optional SHA-256.
`table_info` stores table name, semantic role, and row count.
`column_mapping` always exists and is empty when source headers required no
mapping. `validation_issue` always exists in metadata v1 and records
non-fatal validation findings; its row count must match metadata.

The namespace contains no biological facts. `information_schema` is not used
because it is the SQL system catalog.

## Lazy and atomic write boundary

DuckDB publication streams each lazy relation through a short-lived staging
Parquet, loads it column-wise, writes metadata, checks the database, closes the
connection, and atomically replaces the destination. Transfer Parquet files
are implementation details and never publication artifacts.

`if_exists="fail"` is the default. `"replace"` preserves the previous target
until a complete new staging file is ready.

## Publication boundary

`TidyDataset` has no directory writer and resources expose no `write_tidy()`
compatibility surface. Publications always use an explicit single-file
`write_duckdb(file_out)` destination, where `file_out` denotes the positional
destination value. Metadata schema v1 is the first and only
supported contract. It may include
`bioextract.source_schema_version` only when the upstream source declares an
authoritative schema label. `bioextract.release_version` remains optional.
Writers never infer either value from paths, basenames, archives, timestamps,
or modification times. Legacy development metadata shapes are rejected; v1
never falls back to or dual-writes `bioextract.schema_version`.
An internal `build_id_prefix` may include a path stem for a human-readable
execution label, but that label is never release or source-schema evidence.
