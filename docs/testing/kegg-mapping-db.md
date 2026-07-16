# KEGG Mapping DB Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

This standard covers local KEGG organism mapping inputs, identifier selection,
many-to-many expansion, and tidy publication. The current contract is defined
in [KEGG Mapping DB Architecture](../architecture/kegg-mapping-db.md).

## Automated Tests

[`tests/test_kegg_mapping.py`](../../tests/test_kegg_mapping.py) verifies:

- normalization of UniProt, KEGG gene, KO, pathway, and optional NCBI Gene ID
  inputs
- many-to-many mapping expansion and stable output columns
- single and grouped selection with unmapped reporting
- nullable columns for optional source files
- flat `mapping.parquet` and manifest output
- input-ID kind, snapshot kind, and organism-code validation

The tests use compact local TSV fixtures. They do not call the KEGG API, infer
missing mappings, or calculate enrichment statistics.

## Verification

Run the focused contract with:

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/test_kegg_mapping.py
```

For a real organism snapshot, additionally compare observed source namespaces,
row counts, and a deterministic identifier sample with the raw files before
publishing `mapping.parquet`.
