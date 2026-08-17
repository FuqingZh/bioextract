# ReactomeDatabase Test Standard

Version: v1.2
Date: 2026-08-18
Status: current

## Scope

The suite verifies official-format parsing, the twelve mapping roles, explicit
pathway/reaction dimensions, namespace normalization, human Complex/EWAS
relations, GMT archive safety, species behavior, metadata-v2 publication, and
DuckDB reopen parity. It does not test Reactome web services or enrichment
statistics.

## Fixtures

Mapping fixtures are separate, headerless six-field TSV files. They include an
exact duplicate, evidence-distinct rows, literal quotes, an NCBI non-decimal
identifier, ChEBI prefixed/decimal inputs, and GtoP numeric inputs.

Complex and EWAS fixtures use their exact ordered headers:

```text
complex  pathway  top_level_pathway
ewas     pathway  top_level_pathway
```

They cover an entity prefix that is not `HSA`, exact duplicate removal, missing
metadata endpoints, and an invalid top-level ancestor. GMT fixtures are zip
archives with one `ReactomePathways.gmt` member and cover variable-width rows,
duplicate memberships, opaque symbols, extra members, malformed UTF-8, and
conflicting labels.

All tests use `tmp_path`; they do not write to CephFS or the formal resource
tree.

## Unit And Integration Contract

- Headerless mapping parsing disables quoting, rejects ragged/empty records,
  deduplicates only exact six-field rows, and preserves evidence differences.
- Headered entity parsing requires exact header order and records semantic
  source-to-public column lineage.
- GMT parsing rejects unsafe archive shapes, streams the expected member,
  preserves labels/symbols, and rejects one pathway with multiple labels.
- `from_files()` accepts partial role combinations; absent roles fail at their
  capability boundary rather than being synthesized.
- Pathway `pathway_level=None` resolves to lowest-level; reaction selections
  reject any pathway level. Invalid ChEBI/GtoP identifiers fail before lookup.
- NCBI non-decimal official identifiers remain selectable.
- Existing UniProt selected output remains byte-compatible; reaction output is
  target-specific and source-column-specific.
- Human entity/GMT relations return all rows unscoped or for Homo sapiens and
  an empty stable-schema relation for other species.
- Mapping closure, entity top-level ancestry, cycles, and metadata endpoint
  warnings are conditional on the supporting source roles.

## Publication And Reopen

Partial and complete fixtures verify that `build_tidy()` and `write_duckdb()`
publish exactly the available biological roles plus the five metadata-v2
relations. The complete v0.5 identity is:

```text
reactome-mapping-files-v5
reactome-mapping-v0.5
```

Tests also verify:

- exact role/table/media inventories, including `application/zip` for GMT;
- exact ordered physical schemas and six expected entity column mappings;
- visible warnings without dropping source rows;
- source/reopened whole-resource, species, selection, grouped, unmatched,
  entity, and GMT parity;
- no old `protein_pathway` relation or v0.4-and-earlier reader;
- non-negative recorded table counts are trusted without biological recount;
- read-only independent connections and replacement invalidation.

## Bounded v96 Smoke

External snapshot smoke is opt-in and uses the already resolved concrete raw
subtree. It never recursively scans `/cephfs_data` and writes only to a unique
temporary directory.

The complete v0.5 smoke records all 17 biological tables, source row counts,
exact-duplicate counts, source/publication semantic hashes, four namespace
closure results, selected/grouped probes, human/nonhuman scope probes, GMT
member shape, validation warnings, output hash, elapsed time, and peak RSS only
for operations that actually completed. Fixture results do not substitute for
skipped external evidence.

## Commands

```console
BIOEXTRACT_TEST_THREADS=1 pdm run pytest \
  tests/unit/reactome \
  tests/integration/reactome \
  tests/contract/resources/reactome \
  tests/contract/api/test_signatures.py -q

BIOEXTRACT_TEST_THREADS=1 pdm run check
```

The complete gate is required for this public API and publication-schema
change. Formal replacement, package release, and catalog activation are
separate release operations.
