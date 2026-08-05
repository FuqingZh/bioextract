# KEGG Mapping DB Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

This standard covers local KEGG organism mapping inputs, identifier selection,
many-to-many expansion, and tidy publication. The current contract is defined
in [KEGG Mapping DB Architecture](../architecture/kegg-mapping-db.md).

## Automated Tests

[`tests/integration/kegg/test_mapping.py`](../../tests/integration/kegg/test_mapping.py)
verifies:

- normalization of UniProt, KEGG gene, KO, pathway, and optional NCBI Gene ID
  inputs
- many-to-many mapping expansion and stable output columns
- stable `snake_case` mapping fields and an empty column-provenance inventory
  because the raw mapping files are headerless and the fields are derived
- single and grouped selection with unmapped reporting
- nullable columns for optional source files
- one-table `mapping` DuckDB publication and embedded source inventory
- metadata/profile/schema/table-role validation, read-only connection freshness,
  atomic `if_exists`, and source-to-reopened selection parity
- input-ID kind, snapshot kind, and organism-code validation

The tests use compact local TSV fixtures. They do not call the KEGG API, infer
missing mappings, or calculate enrichment statistics.

## Verification

Run the focused contract with:

```bash
pdm run pytest tests/integration/kegg/test_mapping.py \
  tests/contract/resources/kegg/test_mapping_brite_publication_contract.py \
  tests/contract/api
```

For a real organism snapshot, additionally compare observed source namespaces,
row counts, and a deterministic identifier sample with the raw files before
publishing `tidy/data.duckdb`, then reopen it and repeat the selection through
the publication-backed handle.
