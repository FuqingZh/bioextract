# UniProtDatabase Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.uniprot.UniProtDatabase` provides path-first access to UniProt
`idmapping_selected` resources. The resource is large, so the database handle
must not load data during construction. It supports raw UniProt selected
mapping files, normalized parquet files, and hive parquet datasets from
external or legacy sources.

The MVP covers:

- raw `idmapping_selected.tab` and `idmapping_selected.tab.gz`
- normalized single `mapping.parquet`
- hive parquet dataset directories partitioned by `TaxId`
- `with_taxids(*taxids)` scoped extraction
- full normalized mapping extraction
- single parquet tidy writing
- UniProt `.dat(.gz)` flat-file parsing for eggNOG xref extraction
- UniProt `.dat(.gz)` flat-file parsing for curated Swiss-Prot subcellular
  location comments

It intentionally does not cover:

- online UniProt ID mapping service calls
- broad per-crossref public extraction APIs beyond the explicitly supported
  eggNOG and subcellular-location paths
- enrichment statistics
- GO cellular-component inference from UniProtKB comment text

## Raw Columns

The local `idmapping_selected.tab.gz` has 22 tab-separated columns without a
header:

```text
UniProtId
UniProtEntryName
GeneId
RefSeq
GI
PDB
GO
UniRef100
UniRef90
UniRef50
UniParc
PIR
TaxId
MIM
UniGene
PubMed
EMBL
EMBLCDS
Ensembl
EnsemblTranscript
EnsemblProtein
AdditionalPubMed
```

All columns are normalized as strings.

## Public API

```python
from bioextract.uniprot import UniProtDatabase

db = UniProtDatabase.from_files(
    id_mapping="idmapping_selected.tab.gz",
)

df_hsa = db.with_taxids("9606").extract_mapping()

db.with_taxids("9606", "10090").write_parquet("out/uniprot.parquet")
db.write_parquet("out/uniprot-all.parquet", allow_all_taxa=True)
```

The same constructor accepts tidy outputs:

```python
db = UniProtDatabase.from_files(id_mapping="out/uniprot")
df_hsa = db.with_taxids("9606").extract_mapping()
```

## Construction Checks

`from_files()` performs only lightweight checks:

- path exists
- path type is supported
- file inputs pass configured size limits
- hive dataset directories contain at least one parquet file

It does not read the raw data or collect parquet schemas. Schema validation is
done by `validate_schema()`, `extract_mapping()`, and `write_parquet()`.

## Publication

`write_parquet(path)` publishes one canonical idmapping relation with
footer provenance. It requires `allow_all_taxa=True` when no taxids are selected,
because all-taxa export may scan the entire 9 GB raw gzip file.
`if_exists="fail"` protects an existing file; `"replace"` publishes through a
staging file.

For UniProt knowledge-base flat files, the implemented helper path is:

```python
db = UniProtDatabase.from_dat(
    path="uniprot_sprot.dat.gz",
    source_database="Swiss-Prot",
)
result = db.write_eggnog_xref_parquet("out/uniprot_eggnog_xref.parquet")
```

That path emits a canonical `mapping.parquet` with:

```text
UniProtId
PrimaryUniProtId
IsPrimaryAccession
EggnogOgId
EggnogLevel
SourceDb
```

Swiss-Prot subcellular location comments use the same flat-file constructor:

```python
db = UniProtDatabase.from_dat(
    path="uniprot_sprot.dat.gz",
    source_database="sprot",
)

df_subcell = db.extract_subcellular_location()
result = db.write_subcellular_location_parquet(
    "out/uniprot_subcellular_location.parquet"
)
```

That path emits one Parquet with one row per
`UniProt accession x subcellular location text x evidence`:

```text
UniProtId
PrimaryUniProtId
UniProtEntryName
GeneName
ProteinName
SubcellularLocation
SubcellularLocationNote
EvidenceCode
EvidenceSource
EvidenceId
SourceDb
```

The extractor is deliberately conservative. It preserves curated UniProtKB
`CC   -!- SUBCELLULAR LOCATION:` annotation text and ECO evidence, but it does
not map locations to GO terms or interpret missing comments as negative
localization evidence.

## Implementation Notes

Raw TSV and parquet inputs are scanned lazily. Tidy writing uses Polars
`sink_parquet()`, so large writes do not need to collect the full table in
memory before writing. Hive parquet dataset reading remains supported for
compatibility, but `write_parquet()` publishes one selected relation rather than `TaxId=` partitioned output
because UniProt all-taxa data has very high `TaxId` cardinality.

The shared tidy contract has also changed since the first draft:

- report assets are dataclasses
- manifest asset `sha256` is optional
- `row_count` is no longer part of manifest metadata
- hashing is opt-in through `should_hash_assets=True`
