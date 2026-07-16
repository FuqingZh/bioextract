# GO Term Selection Plan

Version: v1.0
Date: 2026-06-15
Status: superseded by [GO and KEGG Tidy Architecture](../architecture/go-kegg-tidy.md)

## Goal

Add a term-selection layer to `GoDb` that works for both full GO OBO snapshots
and GO subset OBO snapshots such as `goslim_generic.obo`.

The design treats GO slim as OBO subset data, not as a separate database type.
`GoDb.from_obo()` remains the single constructor for GO OBO files. Selection
methods expose stable views over the parsed snapshot.

## Scope

This plan covers:

- preserving OBO subset facts in tidy output
- selecting terms by GO ID, namespace, and subset membership
- listing subset definitions available in the current snapshot
- keeping `extract_subcell()` as a compatibility convenience over term
  selection

This plan does not cover:

- GO GAF, GPAD, or protein/gene annotation selection
- `select_ids()` or `select_groups()` for biological entity IDs
- enrichment statistics
- projecting arbitrary GO annotations to slim ancestors
- full-GO and subset-GO cross-snapshot validation

Those behaviors require an annotation or graph-projection layer separate from
basic OBO term selection.

## Data Model

Extend GO OBO parsing so `[Term]` records retain `subset:` values.

Add these tidy frames:

```text
subset_membership.parquet
go_id
subset_id
```

```text
subset_definition.parquet
subset_id
subset_name
```

Optionally add later, if needed for fuller OBO fidelity:

```text
typedef.parquet
typedef_id
typedef_name
namespace
xref_text
parent_typedef_id
```

Header metadata such as `data-version`, `ontology`, and
`owl:versionInfo` should be stored in source metadata or manifest content
rather than forced into the term table.

The GO tidy schema version must be bumped when these frames are added.

## Public API

Add `GoDb.select_terms()` as the main user-facing query:

```python
from typing import Literal

GoNamespace = Literal[
    "biological_process",
    "cellular_component",
    "molecular_function",
]

def select_terms(
    self,
    *,
    term_ids: Iterable[str] | None = None,
    namespace: GoNamespace | None = None,
    subset_id: str | GoSubsetId | None = None,
    include_obsolete: bool = False,
    should_resolve_alt_ids: bool = True,
) -> pl.DataFrame:
    ...
```

The public `namespace` type should be a `Literal` so common IDEs expose the
three valid values directly. Internally, normalize it to a `StrEnum` or an
equivalent canonical value before filtering.

`subset_id` remains open-ended because OBO subset IDs are data values, not a
closed interface. Provide common constants or a `StrEnum` for discoverability:

```python
db.select_terms(subset_id=GoSubsetId.GOSLIM_GENERIC)
db.select_terms(subset_id="goslim_generic")
```

Add `GoDb.list_subsets()`:

```python
def list_subsets(self) -> pl.DataFrame:
    ...
```

Return columns:

```text
subset_id
subset_name
num_terms
```

This gives users a data-driven way to discover subset IDs present in the
current OBO snapshot.

## Selection Semantics

`select_terms()` returns a stable term view. Conditions are combined with AND:

- `term_ids` filters to requested GO IDs
- `namespace` filters the term namespace
- `subset_id` filters to terms with matching subset membership
- `include_obsolete=False` excludes obsolete terms

Default output is one row per selected primary GO term:

```text
go_id
term_name
namespace
definition
is_obsolete
comment
```

If `subset_id` is provided, include:

```text
subset_id
```

Do not expand all subset memberships by default. Full subset membership remains
available through `build_tidy().frames["subset_membership"]`.

When `should_resolve_alt_ids=True`, input alternate GO IDs should resolve to
primary GO IDs through `alt_id`. Include `input_go_id` when `term_ids` is
provided so callers can see canonicalization.

Missing `subset_id` values should return an empty DataFrame rather than raise.
Strict validation can be added later if a caller needs it.

## Relationship To Existing API

Keep `build_tidy()` as the complete fact-layer API. It should parse the current
OBO snapshot and expose all tidy frames, including subset frames.

Keep `extract_subcell()` for compatibility, but implement it conceptually as a
convenience view over:

```python
select_terms(namespace="cellular_component")
```

It may rename `term_name` to `subcell_name`, but `select_terms()` itself should
keep ontology-level names.

Do not add `kind_db` to `GoDb.from_obo()`. Full GO OBO and GO subset OBO files
share the same format. Their differences belong in metadata and query
parameters, not in constructor branching.

Do not add `select_ids()` or `select_groups()` to OBO-only `GoDb`. In this
repository those names should remain reserved for selection by biological input
IDs against mapping resources. GO can grow those methods later only after an
annotation source such as GAF, GPAD, or another gene/protein-to-GO mapping is
introduced.

## Implementation Steps

1. Extend parser models with term subsets and header subset definitions.
2. Add column buffers and schemas for `subset_membership` and
   `subset_definition`.
3. Add those frames to `build_tidy_frames()` and `ASSET_SPECS`.
4. Bump the GO tidy schema version.
5. Add internal namespace normalization from public `Literal` values to the
   canonical namespace representation.
6. Add `GoSubsetId` common constants for high-use subset IDs such as
   `goslim_generic`, while still accepting arbitrary strings.
7. Add `GoDb.select_terms()` with cached tidy-frame reuse.
8. Add `GoDb.list_subsets()`.
9. Rework `extract_subcell()` to call or mirror `select_terms()` behavior.
10. Update README and GO architecture docs after tests define the final output
    contract.

## Test Contract

Add focused tests covering:

- `goslim_generic.obo`-style subset OBO parses `subset_membership`
- `select_terms(subset_id="goslim_generic")` returns one row per term
- `namespace` and `subset_id` combine with AND semantics
- `list_subsets()` returns subset definitions and term counts
- terms with multiple subset memberships still appear once in `select_terms()`
- obsolete terms are excluded by default and included when requested
- alternate GO IDs resolve to primary GO IDs when `should_resolve_alt_ids=True`
- `extract_subcell()` remains compatible with existing output expectations

Use a compact synthetic OBO fixture for unit tests, then validate against the
real `goslim_generic.obo` snapshot when local resource access is available.
