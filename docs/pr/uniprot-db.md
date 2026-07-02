# PR Plan: UniprotDb Integration

## Summary

- Add `bioextract.uniprot.UniprotDb`.
- Support raw UniProt `idmapping_selected.tab(.gz)`.
- Support normalized single parquet and hive parquet dataset inputs.
- Support taxid-scoped extraction through `with_taxids(*taxids)`.
- Write tidy outputs as canonical single parquet files.
- Require `should_allow_all=True` for all-taxid writing.
- Use `policy_existing="error" | "overwrite" | "skip"` for non-empty outputs.
- Support `level_compression` for zstd parquet compression tuning.

## Files

Add:

```text
src/bioextract/uniprot/__init__.py
src/bioextract/uniprot/constant.py
src/bioextract/uniprot/uniprot.py
src/bioextract/uniprot/util.py
tests/test_uniprot.py
docs/architecture/uniprot-db.md
docs/pr/uniprot-db.md
docs/testing/uniprot-db.md
```

Update:

```text
src/bioextract/__init__.py
README.md
```

## Public Usage

```python
from bioextract.uniprot import UniprotDb

db = UniprotDb.from_files(
    file_idmapping_selected="idmapping_selected.tab.gz",
)

df_hsa = db.with_taxids("9606").extract_mapping()
db.with_taxids("9606", "10090").write_tidy("out/uniprot")
db.write_tidy("out/uniprot-all", should_allow_all=True)
```

Read tidy output:

```python
db = UniprotDb.from_files(file_idmapping_selected="out/uniprot")
df_hsa = db.with_taxids("9606").extract_mapping()
```

## Compatibility

- Additive public API only.
- No change to existing resource DBs.
- No new runtime dependency beyond Polars.
- Tidy output uses `mapping.parquet`; hive parquet dataset reading remains
  available for external or legacy inputs.

## Tests

Run:

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/test_uniprot.py
pdm run precommit
```

Heavy real-data all-taxid writes should not be part of unit tests.
