# KEGG Mapping Db Architecture

## Goal

`bioextract.kegg.KeggDb` should support local KEGG organism mapping snapshots
in addition to the existing BRITE JSON tidy path. The mapping path turns
explicit KEGG raw files into a stable annotation table for downstream
proteomics enrichment.

The first version covers:

- explicit file-based construction, with no hidden directory naming contract
- UniProt, NCBI GeneID, and KEGG gene input selection
- single and grouped selections
- one wide `mapping.parquet` tidy artifact
- optional manifest writing through the shared `TidyDataset` contract

The first version intentionally does not cover:

- full-organism batch publishing
- online KEGG API calls
- enrichment p-value calculation
- pathway hierarchy parsing beyond map-ID derivation
- KO-to-reference-pathway inference when `gene_pathway.tsv` is absent

## Raw Inputs

`KeggDb.from_mapping_files()` accepts exact file paths:

```text
conv_uniprot.tsv
gene_ko.tsv
gene_pathway.tsv
gene_list.tsv
conv_ncbi_geneid.tsv
```

`conv_uniprot.tsv`, `gene_ko.tsv`, and `gene_pathway.tsv` are required for the
first version because the proteomics path needs UniProt-to-pathway and
UniProt-to-KO mappings. `gene_list.tsv` and `conv_ncbi_geneid.tsv` are optional
metadata and alternate-ID inputs.

The organism code is explicit:

```python
KeggDb.from_mapping_files(..., organism_code="hsa")
```

It is also used to validate that observed KEGG gene IDs belong to the expected
organism namespace when the raw files are non-empty.

## Public API

```python
from bioextract.kegg import KeggDb

db = KeggDb.from_mapping_files(
    file_conv_uniprot="conv_uniprot.tsv",
    file_gene_ko="gene_ko.tsv",
    file_gene_pathway="gene_pathway.tsv",
    organism_code="hsa",
    file_gene_list="gene_list.tsv",
    file_conv_ncbi_geneid="conv_ncbi_geneid.tsv",
)

df_mapping = db.extract_mapping()

selection = db.select_ids(["P12345", "Q9Y243"], kind_input_id="uniprot")
df_selected = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()

grouped = db.select_groups(
    {"up": ["P12345"], "down": ["Q9Y243"]},
    kind_input_id="uniprot",
)
df_grouped = grouped.extract_mapping()
```

The accepted `kind_input_id` values are:

```text
uniprot
ncbi_geneid
kegg_gene
```

The name follows the global `structure_role` convention: `kind_...` identifies
the semantic kind of the input ID values.

## Output Contract

`extract_mapping()` and `mapping.parquet` expose one wide table:

```text
OrganismCode
KeggGeneId
UniProtId
NcbiGeneId
KoId
KeggPathwayId
PathwayMapId
GeneSymbol
GeneDescription
```

`KeggGeneId` keeps the KEGG-native namespace such as `hsa:10458`. Other raw
database prefixes are normalized away:

```text
up:P12345        -> UniProtId=P12345
ncbi-geneid:1    -> NcbiGeneId=1
ko:K00001        -> KoId=K00001
path:hsa00010    -> KeggPathwayId=hsa00010
```

`PathwayMapId` is derived from organism-specific pathway IDs:

```text
hsa00010 -> map00010
```

Many-to-many relationships are represented as multiple rows. This keeps the
artifact simple for downstream joins and lets enrichment callers derive
pathway or KO `term2gene` tables with projections.

Single selections prepend:

```text
InputId
KindInputId
```

Grouped selections prepend:

```text
GroupId
InputId
KindInputId
```

`extract_unmapped_input_ids()` returns `InputId` for single selections and
`GroupId, InputId` for grouped selections.

## Tidy Dataset

`build_tidy()` emits:

```text
mapping.parquet
```

Suggested schema version:

```text
kegg-mapping-v0.1
```

The existing BRITE JSON path keeps its current `pathway.parquet` contract.
`KeggDb` should internally distinguish snapshot kinds and dispatch
`build_tidy()` accordingly.

## Implementation Notes

- Keep the public entrypoint on `KeggDb`; do not add another top-level DB type.
- Keep file paths explicit in `from_mapping_files()`.
- Keep selection output flat, matching Reactome, WikiPathways, StringDB, and
  OmniPath grouped selection behavior.
- Cache the full mapping frame on the `KeggDb` instance after first materialize.
- Use Polars joins and projections rather than manual row loops.
- Raise targeted `ValueError` when mapping-only methods are called on a BRITE
  snapshot.

## Tests

Add focused tests for:

- `from_mapping_files()` with required and optional files
- `extract_mapping()` normalization and many-to-many expansion
- `select_ids(..., kind_input_id="uniprot")`
- `select_ids(..., kind_input_id="ncbi_geneid")`
- `select_ids(..., kind_input_id="kegg_gene")`
- `select_groups(..., kind_input_id="uniprot")`
- unmapped IDs for single and grouped selections
- missing optional files with null output columns
- `build_tidy()` and `write_tidy()` writing `mapping.parquet`
- invalid `kind_input_id`
- mapping APIs called on a BRITE snapshot
