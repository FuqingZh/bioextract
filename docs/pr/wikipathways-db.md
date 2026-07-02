# PR Plan: WikiPathwaysDb Integration

## Summary

- Add `bioextract.wikipathways.WikiPathwaysDb`.
- Parse local WikiPathways GMT files.
- Emit `pathway`, `term2gene`, and `term2name` frames.
- Support single and grouped Entrez Gene ID selections.
- Support unmapped input ID reporting.
- Support flat tidy parquet writing with optional manifest.

## Data Contract

Raw input:

```text
wikipathways-YYYYMMDD-gmt-{Species}.gmt
```

GMT row contract:

```text
PathwayName%Collection%WikiPathwaysId%Species<TAB>Url<TAB>EntrezGeneId...
```

The gene IDs are treated as NCBI Entrez Gene IDs. No ID conversion is attempted
in v0.1.

## Files

Add:

```text
src/bioextract/wikipathways/__init__.py
src/bioextract/wikipathways/constant.py
src/bioextract/wikipathways/util.py
src/bioextract/wikipathways/wikipathways.py
tests/test_wikipathways.py
docs/architecture/wikipathways-db.md
docs/pr/wikipathways-db.md
docs/testing/wikipathways-db.md
```

Update:

```text
src/bioextract/__init__.py
README.md
```

## Public Usage

```python
from bioextract.wikipathways import WikiPathwaysDb

db = WikiPathwaysDb.from_gmt(
    "wikipathways-20260510-gmt-Homo_sapiens.gmt",
    species="Homo sapiens",
)

df_term2gene = db.extract_term2gene()
df_term2name = db.extract_term2name()
df_mapping = db.select_ids(["2687", "435"]).extract_mapping()
db.write_tidy("out/wikipathways-hsa", should_write_manifest=True)
```

## Compatibility

- Additive public API only.
- No change to GO, KEGG, STRINGdb, OmniPath, or Reactome behavior.
- No new runtime dependency beyond Polars.
- Output files are flat parquet assets.

## Tests

Run:

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/test_wikipathways.py
PYTHONPATH=src .venv/bin/python -m pytest
```

Real-data smoke:

```bash
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path
from bioextract.wikipathways import WikiPathwaysDb

file_gmt = Path("/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/wikipathways/gmt/2026-05-10/raw/Homo_sapiens/wikipathways-20260510-gmt-Homo_sapiens.gmt")
db = WikiPathwaysDb.from_gmt(file_gmt, species="Homo sapiens")
print(db.extract_pathway().shape)
print(db.extract_term2gene().shape)
print(db.extract_term2name().shape)
db.write_tidy("/tmp/bioextract-wikipathways-hsa", should_write_manifest=True)
PY
```
