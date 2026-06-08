# Tidy Dataset Contract

## Goal

`bioextract._tidy` defines the shared write contract for resource snapshots
that emit parquet artifacts plus an optional manifest.

This document is the authority for:

- `TidyDataset`
- `TidyAsset`
- `TidyReportAsset`
- `TidyManifestAsset`
- `TidyWriteReport`
- `manifest.json` asset metadata shape

It does not define any resource-specific biological schema.

## Current Contract

`TidyDataset` is a lazy-write boundary:

- `frames` must be `pl.LazyFrame`
- `write()` persists frames with `sink_parquet()`
- manifest writing is optional
- asset hashing is optional

The implemented public shape is:

```python
report = dataset.write(
    dir_out,
    should_write_manifest=True,
    should_hash_assets=False,
)
```

## Asset Types

Report assets are dataclasses:

```python
TidyReportAsset(
    path="mapping.parquet",
    kind="canonical",
    is_optional=False,
)
```

Manifest assets are also dataclasses:

```python
TidyManifestAsset(
    path="mapping.parquet",
    kind="canonical",
    sha256=None,
    is_optional=False,
)
```

`TidyWriteReport.assets` returns `tuple[TidyReportAsset, ...]`.

`manifest.json` remains JSON-like and serializable. Manifest asset dataclasses
are converted with `asdict()` at the write boundary.

## Output Contract

All tidy writers emit flat outputs by default:

```text
out/
  <asset>.parquet
  manifest.json
```

There is no default `canonical/` or `derived/` subdirectory split.

The manifest contains:

- `build_id`
- `schema_version`
- `generated_at`
- `sources`
- `assets`

Each manifest asset contains:

- `path`
- `kind`
- `sha256`
- `is_optional`

`row_count` is intentionally no longer part of the manifest contract.

## Hashing Policy

`should_hash_assets=False` is the default for large-resource practicality.

When hashing is disabled:

- manifest asset `sha256` is `null`
- source file metadata still records `path`, `bytes`, and `media_type`

When hashing is enabled:

- each written parquet asset receives a SHA256 digest

The default avoids a second full-file read for large artifacts such as
UniProt, InterPro, and eggNOG outputs.

## Lazy Boundary Rule

The shared direction is:

- construction should stay path-first
- raw table scanning should stay lazy whenever the source format allows it
- `write_tidy()` should not materialize full DataFrames only to write parquet

Small non-tabular parsers remain allowed to do eager parse work when the source
format does not map cleanly to `scan_*()` APIs. Examples include:

- GO OBO stanza parsing
- WikiPathways GMT header parsing
- UniProt `.dat` flat-file record parsing

For those formats, the expectation is still:

- keep eager parsing tightly scoped to the format boundary
- convert to lazy frames as early as practical after parsing

## Resource-Specific Exceptions

Some resources keep custom `write_tidy()` implementations because they need
staging or streaming policy beyond the generic writer:

- `UniprotDb.write_tidy()` stages output before publish and monitors resources
- `InterProDb.write_tidy()` writes a single canonical parquet directly from a
  lazy join plan
- `EggnogDb.write_tidy()` expands SQLite into a temporary TSV before lazy scan

These implementations should still conform to the same external tidy contract.
