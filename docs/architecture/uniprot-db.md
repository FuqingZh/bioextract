# UniprotDb Architecture

## Goal

`bioextract.uniprot.UniprotDb` provides path-first access to UniProt
`idmapping_selected` resources. The resource is large, so the database handle
must not load data during construction. It should support raw UniProt selected
mapping files and the hive parquet datasets written by bioextract.

The MVP covers:

- raw `idmapping_selected.tab` and `idmapping_selected.tab.gz`
- normalized single `mapping.parquet`
- hive parquet dataset directories partitioned by `TaxId`
- `with_taxids(*taxids)` scoped extraction
- full normalized mapping extraction
- hive parquet tidy writing by default

It intentionally does not cover:

- UniProt `.dat` parsing
- online UniProt ID mapping service calls
- per-crossref public extraction APIs
- enrichment statistics

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
from bioextract.uniprot import UniprotDb

db = UniprotDb.from_files(
    file_idmapping_selected="idmapping_selected.tab.gz",
)

df_hsa = db.with_taxids("9606").extract_mapping()

db.with_taxids("9606", "10090").write_tidy("out/uniprot")
db.write_tidy("out/uniprot-all", should_allow_all=True)
```

The same constructor accepts tidy outputs:

```python
db = UniprotDb.from_files(file_idmapping_selected="out/uniprot")
df_hsa = db.with_taxids("9606").extract_mapping()
```

## Construction Checks

`from_files()` performs only lightweight checks:

- path exists
- path type is supported
- file inputs pass configured size limits
- hive dataset directories contain at least one parquet file

It does not read the raw data or collect parquet schemas. Schema validation is
done by `validate_schema()`, `extract_mapping()`, and `write_tidy()`.

## Tidy Output

Hive parquet writing is the default:

```text
out/
  TaxId=9606/
    00000000.parquet
  TaxId=10090/
    00000000.parquet
  manifest.json
```

`write_tidy()` requires `should_allow_all=True` when no taxids are selected,
because all-taxa export may scan the entire 9 GB raw gzip file.

Non-hive output is available only for exactly one selected TaxId:

```python
db.with_taxids("9606").write_tidy("out/hsa", should_write_hive=False)
```

which writes:

```text
out/hsa/
  mapping.parquet
```

## Implementation Notes

Raw TSV and parquet inputs are scanned lazily. Hive writing uses Polars
`sink_parquet(PartitionBy(...))`, so multi-taxid and all-taxid writes do not
need to collect the full table in memory before partitioning.
