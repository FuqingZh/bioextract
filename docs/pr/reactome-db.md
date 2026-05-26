# PR Plan: ReactomeDb Integration

## Summary

- Add `bioextract.reactome.ReactomeDb` for local Reactome mapping snapshots.
- Support UniProt-to-pathway annotation extraction.
- Support species-scoped `term2gene`, `term2name`, and pathway relation frames.
- Support single-query and grouped selections with unmapped ID reporting.
- Add optional tidy parquet output through the shared `TidyDataset` writer.
- Allow STRINGdb-style partial file combinations with targeted feature-level
  errors when a required file is missing.
- Document Reactome as a local enrichment-input resource, not as an online
  Reactome Analysis Service wrapper.

## Public Usage

```python
from bioextract.reactome import ReactomeDb

db = ReactomeDb.from_files(
    file_uniprot2reactome="UniProt2Reactome.txt",
    file_pathways="ReactomePathways.txt",
    file_relations="ReactomePathwaysRelation.txt",
)

selection = (
    db.with_species("Homo sapiens")
    .select_ids(["P04637", "Q9Y243", "MISSING"])
)

df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()

df_term2gene = db.with_species("Homo sapiens").extract_term2gene()
df_term2name = db.with_species("Homo sapiens").extract_term2name()
```

Tidy use:

```python
tidy = db.build_tidy()
tidy.write("out/reactome-v96", should_write_manifest=True)
```

## Files

Add:

```text
src/bioextract/reactome/__init__.py
src/bioextract/reactome/constant.py
src/bioextract/reactome/reactome.py
src/bioextract/reactome/util.py
tests/test_reactome.py
```

Update:

```text
src/bioextract/__init__.py
README.md
```

## Data Contract

Supported raw input files:

```text
UniProt2Reactome.txt
ReactomePathways.txt
ReactomePathwaysRelation.txt
```

At least one file is required. Individual capabilities depend only on their
backing files: mapping and `term2gene` use `UniProt2Reactome.txt`; `term2name`
uses `ReactomePathways.txt`; unscoped relations use
`ReactomePathwaysRelation.txt`; species-scoped relations also require
`ReactomePathways.txt`.

Public extraction columns:

```text
InputId
UniProtId
ReactomePathwayId
PathwayName
EvidenceCode
Species
ReactomeUrl
```

Grouped extraction prepends `GroupId`.

Tidy schema version:

```text
reactome-mapping-v0.1
```

Tidy output:

```text
mapping.parquet
pathway.parquet
relation.parquet
term2gene.parquet
term2name.parquet
```

## Compatibility

- Additive public API only.
- No change to STRINGdb, OmniPath, GO, or KEGG behavior.
- No new runtime dependency beyond Polars.
- Raw Reactome display names are preserved.
- ID selection is UniProt-only in the first version.

## Implementation Steps

1. Add Reactome constants, column names, and required file metadata.
2. Add `ReactomeDb.from_files()` with optional files, path validation, and file
   size limits.
3. Add lazy Polars readers for mapping, pathway, and relation frames.
4. Add `with_species()`, `select_ids()`, and `select_groups()`.
5. Add mapping and unmapped extraction.
6. Add `term2gene`, `term2name`, and species-limited relation extraction.
7. Add `build_tidy()` and `write_tidy()` using `TidyDataset`, emitting only
   derivable frames for partial snapshots.
8. Add README examples.
9. Add unit tests and a small real-data smoke command to the PR notes.

## Open Decisions

- Whether species matching should remain exact or also support case-insensitive
  aliases. Default recommendation: exact matching for v0.1.
- Whether relations should include parent pathways outside the selected species
  metadata. Default recommendation: no; keep both endpoints inside the selected
  species.
- Whether to include `ReactomeUrl` in `term2name`. Default recommendation: no;
  keep URL in mapping and pathway frames only.

## Tests

Run:

```bash
PYTHONPATH=src pytest
```

Focused tests:

```bash
PYTHONPATH=src pytest tests/test_reactome.py
```

Real-data smoke:

```bash
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path
from bioextract.reactome import ReactomeDb

base = Path("/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/reactome/mapping/v96/raw")
db = ReactomeDb.from_files(
    file_uniprot2reactome=base / "UniProt2Reactome.txt",
    file_pathways=base / "ReactomePathways.txt",
    file_relations=base / "ReactomePathwaysRelation.txt",
)

view = db.with_species("Homo sapiens")
print(view.extract_term2gene().shape)
print(view.extract_term2name().shape)
print(view.extract_pathway_relations().shape)
PY
```
