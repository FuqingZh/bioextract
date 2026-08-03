# ReactomeDatabase Test Standard

Version: v1.0
Date: 2026-07-30
Status: current

## Scope

The ReactomeDatabase test suite verifies local raw-file parsing,
selection behavior, species scoping, enrichment input extraction, relation
filtering, tidy writing, and error handling.

It does not test online Reactome services, ReactomePA behavior, or statistical
enrichment calculations.

## Fixtures

Use tiny tab-separated fixtures that preserve the real file shapes.

`UniProt2Reactome.txt`:

```text
P04637	R-HSA-69563	https://reactome.org/PathwayBrowser/#/R-HSA-69563	p53-Dependent G1 DNA Damage Response	Gene Ontology	Homo sapiens
P04637	R-HSA-6798695	https://reactome.org/PathwayBrowser/#/R-HSA-6798695	Neutrophil degranulation	TAS	Homo sapiens
Q9Y243	R-HSA-6798695	https://reactome.org/PathwayBrowser/#/R-HSA-6798695	Neutrophil degranulation	TAS	Homo sapiens
P31749	R-MMU-1257604	https://reactome.org/PathwayBrowser/#/R-MMU-1257604	PIP3 activates AKT signaling	TAS	Mus musculus
```

`ReactomePathways.txt`:

```text
R-HSA-69563	p53-Dependent G1 DNA Damage Response	Homo sapiens
R-HSA-6798695	Neutrophil degranulation	Homo sapiens
R-HSA-1640170	Cell Cycle	Homo sapiens
R-MMU-1257604	PIP3 activates AKT signaling	Mus musculus
```

`ReactomePathwaysRelation.txt`:

```text
R-HSA-1640170	R-HSA-69563
R-HSA-1640170	R-HSA-6798695
R-MMU-000001	R-MMU-1257604
```

## Unit Tests

### Construction

- `from_files()` accepts all three files.
- `from_files()` accepts mapping-only, pathways-only, and relations-only
  snapshots.
- `from_files()` rejects a call with no files.
- `from_files()` rejects missing provided files.
- Missing-file failures happen at the feature boundary. For example,
  `extract_term2name()` without pathways raises a targeted `ValueError`.

### Species Scope

- `with_species("Homo sapiens")` filters mapping and pathway frames.
- `with_species("Mus musculus")` returns mouse-only rows.
- Unknown species returns empty extraction frames, not a crash.
- Species matching trims whitespace.

### Single Selection

- `select_ids([" P04637 ", "", "MISSING"])` trims input IDs and drops blanks.
- `extract_mapping()` returns only selected IDs.
- Duplicate raw mappings are deduplicated by UniProt and pathway ID.
- `extract_unmatched_ids()` reports `MISSING`.
- Empty selection returns empty mapping and unmapped frames with stable columns.

### Grouped Selection

- `select_groups()` preserves group labels.
- Grouped mapping prepends `GroupId`.
- Grouped unmapped output includes `GroupId` and `InputId`.
- The same UniProt ID can appear in multiple groups without being collapsed
  across groups.

### Enrichment Input Frames

- `extract_term2gene()` returns `ReactomePathwayId, UniProtId`.
- `extract_term2gene()` is species-scoped.
- `extract_term2name()` returns `ReactomePathwayId, PathwayName, Species`.
- `extract_term2name()` is deduplicated by pathway ID.

### Relations

- `extract_pathway_relations()` returns parent-child columns.
- Species-scoped relations keep only edges where both endpoints are in
  species-scoped pathway metadata.
- Relations are stable when relation rows reference pathways absent from
  metadata.
- Relations-only snapshots can extract unscoped raw relations.
- Species-scoped relations without pathways raise a targeted `ValueError`.

### Tidy

- `build_tidy().frames` contains `mapping`, `pathway`, `relation`,
  `term2gene`, and `term2name`.
- `write_duckdb(path)` writes canonical relations to one database.
- `_bioextract.table_info` agrees with every physical table row count.
- Partial snapshots write only derivable relations and list only provided
  source files.

## Real-Data Smoke

Use the local v96 snapshot:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/reactome/mapping/v96/raw
```

Expected raw sizes observed during design:

```text
UniProt2Reactome.txt: 322435 rows
ReactomePathways.txt: 23498 rows
ReactomePathwaysRelation.txt: 23612 rows
```

The smoke test should verify:

- files read without missing columns
- Homo sapiens `term2gene` and `term2name` are non-empty
- relation filtering completes
- DuckDB publication succeeds under `/tmp`
- the formal v96 publication has 322,435 `protein_pathway`, 23,498 `pathway`,
  and 23,612 `pathway_relation` rows
- metadata schema v3 has all five `_bioextract` tables and reports
  `validation_status=passed`

Example command:

```bash
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path
from bioextract import ReactomeDatabase

base = Path("/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/reactome/mapping/v96/raw")
db = ReactomeDatabase.from_files(
    uniprot_mapping=base / "UniProt2Reactome.txt",
    pathways=base / "ReactomePathways.txt",
    relations=base / "ReactomePathwaysRelation.txt",
)
view = db.with_species("Homo sapiens")
assert view.extract_term2gene().height > 0
assert view.extract_term2name().height > 0
assert view.extract_pathway_relations().height > 0
view.write_duckdb("/tmp/reactome.duckdb")
PY
```

## Non-Goals

- No live Reactome API calls in unit tests.
- No statistical enrichment p-value tests.
- No ReactomePA output parity test.
- No gene-symbol conversion test until ID mapping is explicitly added.
