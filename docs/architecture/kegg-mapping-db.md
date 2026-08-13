# KEGG Mapping Db Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.kegg.KEGGDatabase` supports local KEGG organism mapping snapshots
in addition to the existing BRITE JSON tidy path. The mapping path turns
explicit KEGG raw files into stable gene, KO, and pathway relationships.

The first version covers:

- explicit file-based construction, with no hidden directory naming contract
- UniProt, NCBI GeneID, and KEGG gene input selection
- single and grouped selections
- one wide `mapping` relation in a metadata-v1 DuckDB publication
- embedded source-role provenance with no sidecar manifest

The first version intentionally does not cover:

- full-organism batch publishing
- online KEGG API calls
- enrichment p-value calculation
- pathway hierarchy parsing beyond map-ID derivation
- KO-to-reference-pathway inference when `gene_pathway.tsv` is absent

## Raw Inputs

`KEGGDatabase.from_mapping_files()` accepts exact file paths:

```text
conv_uniprot.tsv
gene_ko.tsv
gene_pathway.tsv
gene_list.tsv
conv_ncbi_gene.tsv
```

`conv_uniprot.tsv`, `gene_ko.tsv`, and `gene_pathway.tsv` are required for the
first version because together they define UniProt-to-pathway and
UniProt-to-KO mappings. `gene_list.tsv` and `conv_ncbi_gene.tsv` are optional
metadata and alternate-ID inputs.

The organism code is explicit:

```python
KEGGDatabase.from_mapping_files(..., organism_code="hsa")
```

It is also used to validate that observed KEGG gene IDs belong to the expected
organism namespace when the raw files are non-empty.

## Public API

```python
from bioextract import KEGGDatabase

db = KEGGDatabase.from_mapping_files(
    uniprot_conversion="conv_uniprot.tsv",
    gene_ko="gene_ko.tsv",
    gene_pathway="gene_pathway.tsv",
    organism_code="hsa",
    gene_list="gene_list.tsv",
    ncbi_gene_conversion="conv_ncbi_gene.tsv",
)

lf_mapping = db.mappings()

selection = db.select_ids(["P12345", "Q9Y243"], namespace="uniprot")
lf_selected = selection.mappings()
lf_unmapped = selection.unmatched_ids()

grouped = db.select_groups(
    {"up": ["P12345"], "down": ["Q9Y243"]},
    namespace="uniprot",
)
lf_grouped = grouped.mappings()

result = db.write_duckdb("kegg-mapping.duckdb")
published = KEGGDatabase.from_duckdb(result.path)
with published.connect() as connection:
    row_count = connection.sql("SELECT count(*) FROM mapping").fetchone()[0]
```

The accepted `namespace` values are:

```text
uniprot
ncbi_gene
kegg_gene
```

The name follows the global `structure_role` convention: `kind_...` identifies
the semantic kind of the input ID values.

## Output Contract

`mappings()` exposes one wide public lazy relation. Callers choose when to
collect it. These fields are
derived from headerless KEGG inputs, so they are created directly with the
stable public `snake_case` names; no source-header mapping is recorded for
this relation:

```text
organism_code
kegg_gene_id
uniprot_id
ncbi_gene_id
ko_id
kegg_pathway_id
pathway_map_id
gene_symbol
gene_description
```

`kegg_gene_id` keeps the KEGG-native namespace such as `hsa:10458`. Other raw
database prefixes are normalized away:

```text
up:P12345        -> uniprot_id=P12345
ncbi-geneid:1    -> ncbi_gene_id=1
ko:K00001        -> ko_id=K00001
path:hsa00010    -> kegg_pathway_id=hsa00010
```

`pathway_map_id` is derived from organism-specific pathway IDs:

```text
hsa00010 -> map00010
```

Many-to-many relationships are represented as multiple rows. This keeps the
artifact simple for downstream joins and lets enrichment callers derive
pathway or KO `term2gene` tables with projections.

Single selections prepend:

```text
input_id
input_namespace
```

Grouped selections prepend:

```text
group_id
input_id
input_namespace
```

`unmatched_ids()` returns `input_id` for single selections and
`group_id, input_id` for grouped selections.

## Tidy Dataset

`write_duckdb()` publishes one physical relation:

```text
main.mapping
```

The DuckDB relation uses the same `snake_case` columns as
`mappings()`; reopening does not perform an inverse rename and
preserves single/grouped selection and unmatched-ID behavior.

Suggested schema version:

```text
kegg-mapping-v0.1
```

The mapping profile is identified by `kegg-organism-mapping-files-v1`, scope
`mapping`, schema version `kegg-mapping-v0.1`, exact table role `canonical`,
and embedded `bioextract.organism_code`. `from_duckdb()` rejects mismatched
profiles, schemas, roles, table inventories, and row counts. Every `connect()`
call creates a fresh caller-owned read-only connection.

## Implementation Notes

- Keep the public entrypoint on `KEGGDatabase`; do not add another top-level DB type.
- Keep file paths explicit in `from_mapping_files()`.
- Keep selection output flat, matching Reactome, WikiPathways, StringDB, and
  OmniPath grouped selection behavior.
- Keep mapping execution replayable; do not make a full eager mapping cache part
  of the public contract.
- Use Polars joins and projections rather than manual row loops.
- Raise targeted `ValueError` when mapping-only methods are called on a BRITE
  snapshot.

## Tests

Focused tests cover:

- `from_mapping_files()` with required and optional files
- `mappings()` normalization and many-to-many expansion
- `select_ids(..., namespace="uniprot")`
- `select_ids(..., namespace="ncbi_gene")`
- `select_ids(..., namespace="kegg_gene")`
- `select_groups(..., namespace="uniprot")`
- unmapped IDs for single and grouped selections
- missing optional files with null output columns
- `build_tidy()`, `write_duckdb(path)`, and metadata-v1 reopening
- fresh read-only native connections and atomic `if_exists` behavior
- invalid `namespace`
- mapping APIs called on a BRITE snapshot
