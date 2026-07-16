# GO And KEGG Tidy Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

This standard covers GO OBO parsing and selection plus KEGG BRITE JSON parsing
and tidy writing. The current architecture is documented in
[GO and KEGG Tidy Architecture](../architecture/go-kegg-tidy.md).

## Automated Tests

[`tests/test_go_tidy.py`](../../tests/test_go_tidy.py) verifies:

- ontology, relation, synonym, subset, ancestor, and depth frame contracts
- subset discovery and term selection by namespace, subset, primary ID, and
  alternate ID
- obsolete-term defaults and subcellular-component compatibility output
- direct and legacy tidy writers

[`tests/test_kegg_tidy.py`](../../tests/test_kegg_tidy.py) verifies:

- pathway-level and entry/KO parsing variants
- BRITE tidy frame schemas and flat Parquet output
- direct and legacy tidy writers

The tests use compact local fixtures. They do not call remote GO or KEGG
services and do not calculate enrichment statistics.

## Verification

Run the focused contract with:

```bash
PYTHONPATH=src .venv/bin/python -m pytest \
  tests/test_go_tidy.py tests/test_kegg_tidy.py
```

Real-snapshot validation should additionally confirm that produced schemas
match the architecture contract and that all relationship endpoints reference
known terms or pathways under the documented source semantics.
