# UniProtKB Domain Access Plan

Storage status: superseded by the current
[Storage And Publication Convergence Plan](20260803-v1.0-storage-publication-convergence-implementation-plan.md).
The decisions and measurements below remain historical; metadata v3 and the
former single-relation Parquet comparison are not current publication authority.
Any former `from_release` wording is historical guidance, not a current public
constructor contract.

Date: 2026-07-30
Status: implemented

## Accepted boundary

UniProt has two independent products. `idmapping_selected` remains a lazy,
single-relation Parquet product. Reviewed UniProtKB/Swiss-Prot is published as
one relational DuckDB built from exact caller-declared roles:

```python
source = UniProtDatabase.from_knowledgebase(
    entries="uniprot_sprot.dat.gz",
    canonical_sequences="uniprot_sprot.fasta.gz",
    isoform_sequences="uniprot_sprot_varsplic.fasta.gz",
    release_version="2026_01",
)
source.write_duckdb("uniprot.duckdb")
```

`entries` is required. The FASTA roles are optional. Plain and transparently
compressed files are accepted regardless of basename. Directories, archives
containing multiple roles, role discovery, and a `from_release` alias are not
supported.

## Version contract

Four version axes remain separate:

- `bioextract.release_version` is optional official release identity. The
  2026_01 DAT, canonical FASTA, and varsplic FASTA have no global authoritative
  release marker, so this API records it only when supplied by the caller.
- `bioextract.source_schema_profile` is the required bioextract-owned,
  content-validated parser profile. UniProtKB uses
  `uniprotkb-flat-file-v1`.
  Metadata also records independent role profiles for Swiss-Prot DAT,
  canonical FASTA, and varsplic FASTA, including explicit absence.
- `bioextract.source_schema_version` is optional and is absent because these
  inputs do not declare an upstream schema label.
- `bioextract.resource_schema_version` versions the published bioextract
  relations.

File names, parent directories, archive names, timestamps, and gzip metadata
never supply release or schema identity.

Metadata schema v3 introduced `bioextract.resource_schema_version` and required
`bioextract.source_schema_profile`. New v3 publications do not dual-write the
old `bioextract.schema_version`. Readers of legacy metadata v1/v2 require that
old key and treat source schema as unknown.

## Validation and publication

The DAT parser is record-streaming and rejects unterminated records, missing
mandatory `ID`, `AC`, `OX`, or `SQ` facts, duplicate primary accessions,
unreviewed records, and inconsistent ID/SQ/sequence lengths. It validates and
publishes SQ molecular weight and CRC64. Optional canonical FASTA is compared
exactly by identifier and sequence against DAT through a disk-backed index.
Each DAT alternative-product block becomes one ordered entry-context product.
The first IsoId is its main ID; later values are retained as ordered old IDs in
`protein_isoform_identifier`, matching the
[UniProt XML/release-note model](https://www.uniprot.org/release-notes/2019-12-18-release).
Related entries may legitimately repeat the same product in `External` and
owned contexts, as documented by
[UniProt Alternative products](https://www.uniprot.org/help/alternative_products).

Varsplic identifiers resolve through all official IDs to exactly one
`Alternative` owner product; displayed products reuse canonical sequences.
`Sequence=VSP...` publishes ordered, context-keyed
`protein_isoform_variation` rows, normalizes status to `Alternative`, and
requires a materialized sequence for every Alternative product whenever the
varsplic role is supplied.

Temporary TSV relation spools and a SQLite validation index bound Python
memory. The shared writer stages DuckDB, validates metadata/table inventories
and row counts through a reopened read-only connection, and atomically replaces
the destination only after validation.

## Domain access

`from_duckdb()` validates metadata v3, resource identity/schema, relation
inventory, and row-count parity. `connect()` returns a caller-owned native
read-only DuckDB connection.

Supported identifier namespaces are `uniprot`, `entry_name`, `gene_name`,
`gene_id`, `refseq`, `ensembl`, and `isoform_id`. Primary, secondary, main
isoform, and old isoform identifiers may map to multiple current protein
contexts; every valid match is retained. Grouped selection preserves caller
labels. Core extractors cover proteins, accessions, protein and gene names, EC
numbers, GO annotations, cross-references, comments, subcellular locations,
keywords, sequences, isoform products and identifiers, sequence variations,
isoform-to-variation linkage, and unmatched identifiers.

Relations whose semantics are not yet implemented are not published as empty
placeholders. Deferred scope includes full TrEMBL scale acceptance, all feature
types beyond `VAR_SEQ`, full evidence normalization, host-organism projection,
and persistent-index benchmarking.
