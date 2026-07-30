# UniProtDatabase Architecture

Version: v1.0
Date: 2026-07-30
Status: current

`UniProtDatabase` exposes independent idmapping and reviewed UniProtKB
products. Idmapping uses `from_idmapping(path, release_version=None)`,
`scan_mapping()`, `read_mapping()`, and `write_parquet()`. An unscoped eager
read or write requires `allow_all_taxa=True`; lazy scanning remains unscoped by
default.

Swiss-Prot has one raw constructor:

```python
UniProtDatabase.from_knowledgebase(
    entries=...,
    canonical_sequences=None,
    isoform_sequences=None,
    release_version=None,
)
```

Arguments declare exact roles. Content and compression magic determine parser
behavior; basenames and directories do not. The content-validated bundle
profile is `uniprotkb-flat-file-v1`. The files do not expose a global release
or upstream schema label, so release identity is caller-only and
`source_schema_version` is absent.

The streaming publisher uses temporary relation spools plus a disk-backed
sequence/isoform index and commits atomically. DAT is authoritative. Optional
canonical FASTA must match every DAT primary accession and sequence. Optional
varsplic FASTA identifiers must resolve to DAT isoform definitions.

## Published Relations

The `uniprot-knowledgebase-duckdb-v1` schema has exactly these 16 `main`
relations:

| Relation | Purpose |
| --- | --- |
| `protein` | Reviewed entry identity, taxonomy, existence, and version facts |
| `protein_accession` | Ordered primary and secondary accessions |
| `protein_sequence` | Canonical and materialized isoform sequences |
| `protein_isoform` | Ordered entry-context products keyed by their main IsoId |
| `protein_isoform_identifier` | Ordered main and old official IsoIds |
| `protein_sequence_variation` | DAT `VAR_SEQ` features |
| `protein_isoform_variation` | Ordered entry-context product-to-VSP relationships |
| `protein_name` | Ordered recommended, alternative, and submitted names |
| `gene_name` | Ordered official gene-name classes |
| `protein_ec_number` | EC annotations |
| `protein_go_annotation` | GO term, aspect, and evidence annotations |
| `protein_cross_reference` | External database identifiers and isoform scope |
| `protein_comment` | Parsed comment blocks |
| `protein_subcellular_location` | Individual locations and optional notes |
| `protein_keyword` | Ordered keywords |
| `protein_identifier` | Internal namespace index used for selection |

`protein_isoform` is an entry-context product relation keyed by
`primary_accession + isoform_id`; `isoform_id` is the first, main ID in the
official IsoId list. `protein_isoform_identifier` retains every listed ID in
order and marks the main ID. This follows the UniProt XML model, where the first
`id` is the main product ID and later values are old product IDs, while also
preserving the documented case where related entries list the same products in
each entry. [UniProt release 2019_11](https://www.uniprot.org/release-notes/2019-12-18-release)
and [Alternative products help](https://www.uniprot.org/help/alternative_products)
define these semantics.

Official `External` declarations may therefore repeat a product owned and
displayed or materialized by another entry. `Displayed` points to the declaring
entry's canonical sequence. `External` and `Not described` remain unresolved.
Varsplic headers resolve through the identifier relation and update only the
unique `Alternative` owner product.

## Selection And Extractor Schemas

Every matched extractor begins with the stable selection prefix
`GroupId, InputId, InputNamespace, UniProtId`, then adds:

| Extractor | Stable additional columns |
| --- | --- |
| `extract_proteins` | `EntryName, IsReviewed, TaxonId, ProteinExistence, SequenceLength, MolecularWeight, SequenceVersion, EntryVersion` |
| `extract_accessions` | `Accession, AccessionOrder, IsPrimaryAccession` |
| `extract_protein_names` | `NameType, ProteinName, NameOrder` |
| `extract_gene_names` | `NameType, GeneName, NameOrder` |
| `extract_ec_numbers` | `ECNumber` |
| `extract_go_annotations` | `GOId, Aspect, TermName, EvidenceCode, EvidenceSource` |
| `extract_cross_references` | `Database, ExternalId, Properties, IsoformId` |
| `extract_comments` | `CommentId, CommentType, CommentText` |
| `extract_subcellular_locations` | `SubcellularLocation, SubcellularLocationNote` |
| `extract_keywords` | `Keyword, KeywordOrder` |
| `extract_sequences` | `SequenceId, SequenceType, Sequence, Length, CRC64, SHA256` |
| `extract_isoforms` | `IsoformId, IsoformName, IsoformOrder, SequenceStatus, SequenceId` |
| `extract_isoform_identifiers` | `IsoformId, Identifier, IdentifierOrder, IsMain` |
| `extract_sequence_variations` | `VariationId, StartPosition, EndPosition, Note` |
| `extract_isoform_variations` | `IsoformId, VariationId, VariationOrder` |

`extract_unmatched_ids()` instead returns
`GroupId, InputId, InputNamespace, Reason`. Empty selections preserve the
corresponding schema, and every extractor has an explicit domain order after
the stable selection prefix.

`from_duckdb()` validates metadata/resource schemas, exact table inventories,
physical column order/types/nullability, and row counts.
`connect()` is read-only. Selection supports `uniprot`, `entry_name`,
`gene_name`, `gene_id`, `refseq`, `ensembl`, and `isoform_id`, retaining every
canonical match. Extractors expose proteins, accessions, names, EC, GO,
cross-references, comments, locations, keywords, sequences, isoforms, isoform
identifiers, variations, and unmatched inputs.

See the [implementation plan](../implementation-plan/20260730-uniprot-knowledgebase-domain-access.md)
for validation details and deferred scope.
