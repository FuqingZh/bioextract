# PR Notes: GO and KEGG Tidy Integration

## Summary

- Add `bioextract.go.GoDb` for GO OBO tidy datasets.
- Add `bioextract.kegg.KeggDb` for KEGG BRITE JSON tidy datasets.
- Add shared `TidyDataset.write()` parquet writer with optional manifest output.
- Migrate GO ontology and KEGG BRITE tidy internals from `biotidy`.
- Preserve legacy `run_tidy_go_ontology()` and `run_tidy_kegg_brite()` entrypoints
  under the new namespaces.

## Public Usage

```python
from bioextract.go import GoDb
from bioextract.kegg import KeggDb

go_tidy = GoDb.from_obo("go-basic.obo").build_tidy()
go_tidy.write("out/go-basic")
go_tidy.write("out/go-basic-archive", should_write_manifest=True)

kegg_tidy = KeggDb.from_brite_json("br08901.json").build_tidy()
kegg_tidy.write("out/br08901")
```

Convenience wrappers:

```python
GoDb.from_obo("go-basic.obo").write_tidy("out/go-basic")
KeggDb.from_brite_json("br08901.json").write_tidy("out/br08901")
```

## Compatibility

- GO `schema_version`: `go-obo-tidy-v0.1`
- KEGG `schema_version`: `kegg-brite-tidy-v0.1`
- New default parquet output is flat, for example `term.parquet` and
  `pathway.parquet`.
- `manifest.json` is only written when `should_write_manifest=True`.
- No `axiomkit` or `orjson` dependency is required for the migrated code.

## Tests

Run:

```bash
PYTHONPATH=src pytest
```

Focused tests:

```bash
PYTHONPATH=src pytest tests/test_go_tidy.py tests/test_kegg_tidy.py
```
