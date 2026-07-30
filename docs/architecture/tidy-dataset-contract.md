# Materialized Dataset Contract

Version: v1.0
Date: 2026-07-29
Status: current

## Purpose

This contract defines how `bioextract` publishes validated local snapshots. It
is subordinate to the repository-wide
[Domain Access Architecture](20260729-v1.0-domain-access-architecture.md):
materialization is an execution strategy, not the product purpose.

## Output unit

The biological relation structure determines the container:

- one independently useful analytical relation: one Parquet file;
- multiple related relations normally queried together: one DuckDB file;
- an already efficient official query format: direct access unless
  materialization has an identified benefit.

A partial capability of a multi-relation resource retains the same DuckDB
container. Missing inputs produce absent tables, not a different packaging
model and not misleading empty observations.

Canonical writers are:

```python
dataset.write_parquet(path)
dataset.write_duckdb(path, table_names=...)
```

Every writer requires `path`. The library does not infer a semantic
filename.

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

## Parquet provenance

Parquet footer key-value metadata carries:

- `bioextract.metadata_schema_version`;
- `bioextract.resource_name`;
- `bioextract.resource_schema_version`;
- `bioextract.source_schema_profile`;
- optional `bioextract.source_schema_version`;
- optional `bioextract.release_version`;
- optional `bioextract.release_version_source`;
- `bioextract.package_version`;
- `bioextract.generated_at`;
- `bioextract.validation_status`;
- `bioextract.validation_issue_count`;
- `bioextract.sources`;
- optional `bioextract.scope`;
- `bioextract.column_mapping`.

No sidecar file is required. A Parquet publication is complete and
machine-identifiable as one file.

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
mapping. `validation_issue` always exists in metadata v2/v3 and records
non-fatal validation findings; its row count must match metadata.

The namespace contains no biological facts. `information_schema` is not used
because it is the SQL system catalog.

## Lazy and atomic write boundary

Lazy relations remain `pl.LazyFrame` until publication. Parquet output uses
Polars `sink_parquet()`; the public method remains `write_parquet()` because it
also owns validation, provenance, staging, and commit.

DuckDB publication streams each lazy relation through a short-lived staging
Parquet, loads it column-wise, writes metadata, checks the database, closes the
connection, and atomically replaces the destination. Transfer Parquet files
are implementation details and never publication artifacts.

`if_exists="fail"` is the default. `"replace"` preserves the previous target
until a complete new staging file is ready.

## Publication boundary

`TidyDataset` has no directory writer and resources expose no `write_tidy()`
compatibility surface. Publications always use an explicit single-file
`write_parquet(path)` or `write_duckdb(path)` destination.
Metadata schema v3 makes both keys above required. It may include
`bioextract.source_schema_version` only when the upstream source declares an
authoritative schema label. `bioextract.release_version` remains optional.
Writers never infer either value from paths, basenames, archives, timestamps,
or modification times. Legacy metadata v1/v2 uses
`bioextract.schema_version`; v3 never falls back to or dual-writes that name.
An internal `build_id_prefix` may include a path stem for a human-readable
execution label, but that label is never release or source-schema evidence.
