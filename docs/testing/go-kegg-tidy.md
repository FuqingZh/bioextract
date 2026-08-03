# GO And KEGG Tidy Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

This standard covers GO OBO parsing and selection plus KEGG BRITE JSON parsing
and tidy writing. The current architecture is documented in
[GO and KEGG Tidy Architecture](../architecture/go-kegg-tidy.md).

## Automated Tests

[`tests/integration/go/test_ontology.py`](../../tests/integration/go/test_ontology.py)
verifies:

- ontology, relation, synonym, subset, ancestor, and depth frame contracts
- subset discovery and term selection by namespace, subset, primary ID, and
  alternate ID
- obsolete-term defaults and subcellular-component compatibility output
- direct and legacy tidy writers

[`tests/unit/kegg/brite/test_parser.py`](../../tests/unit/kegg/brite/test_parser.py)
and
[`tests/integration/kegg/brite/test_publication.py`](../../tests/integration/kegg/brite/test_publication.py)
verify:

- pathway-level and entry/KO parsing variants
- BRITE tidy frame schemas and one-table DuckDB output
- metadata-v1 profile/table validation and fresh read-only native connections

The tests use compact local fixtures. They do not call remote GO or KEGG
services and do not calculate enrichment statistics.

## Verification

Run the focused contract with:

```bash
pdm run pytest tests/unit/kegg/brite tests/integration/kegg/brite \
  tests/contract/resources/kegg/test_mapping_brite_publication_contract.py
```

Real-snapshot validation should additionally confirm that produced schemas
match the architecture contract and that all relationship endpoints reference
known terms or pathways under the documented source semantics.
