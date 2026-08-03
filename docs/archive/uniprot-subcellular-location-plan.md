# UniProtKB Subcellular Location Plan

Version: v1.0
Date: 2026-06-17
Status: superseded by [UniProtDatabase Architecture](../architecture/uniprot-db.md)

## Goal

Add a local, reproducible extractor for Swiss-Prot curated subcellular location
annotations from UniProtKB flat files such as:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/uniprot/kb/2026_01/raw/knowledgebase/complete/uniprot_sprot.dat.gz
```

The feature should parse `CC   -!- SUBCELLULAR LOCATION:` comments, preserve
ECO evidence where present, retain UniProt entry identity fields, and optionally
join to the existing `idmapping_selected` mapping layer through UniProt
accessions.

This is a UniProtKB comment-annotation extractor. It is not a GO xref parser,
GO CC enrichment implementation, localization predictor, or custom ontology.

## Scientific Boundary

The source annotation is mature and standard: Swiss-Prot subcellular location
comments are curated UniProtKB annotations. The implementation should remain a
conservative extraction layer:

- preserve source text rather than inventing a normalized location ontology
- preserve evidence code and evidence source IDs
- keep accession-level identity explicit
- treat missing Swiss-Prot annotation as unknown coverage, not as absence of
  localization

Do not infer GO CC terms in this extractor. If a later workflow needs GO terms,
add a separate mapping step using an explicit resource such as UniProtKB
Subcellular Location2GO.

## Existing Context

`UniProtDatabase.from_dat()` already accepts UniProtKB `.dat` and `.dat.gz` files.
The current flat-file implementation is narrow: it extracts `DR   eggNOG;`
records into a UniProt-to-eggNOG xref table.

This plan extends the same flat-file handle with subcellular location comment
extraction. It should not change the existing `idmapping_selected` behavior.

## Public API

Add a DataFrame extraction method:

```python
db = UniProtDatabase.from_dat(
    path="uniprot_sprot.dat.gz",
    source_database="sprot",
)

df_subcell = db.extract_subcellular_location()
```

Add a tidy writer:

```python
report = db.write_subcellular_location_tidy(
    "out/subcellular_location",
    should_write_manifest=True,
)
```

The writer emits:

```text
subcellular_location/
  data.parquet
  manifest.json
```

Use `data.parquet` because the directory name carries the semantic asset name.
Do not reuse `mapping.parquet`, which is already associated with idmapping and
xref mapping outputs.

Optional accession mapping should stay explicit:

```python
df_subcell = db_dat.extract_subcellular_location()
df_mapping = db_idmapping.with_taxids("9606").extract_mapping()
df_joined = join_subcellular_location_idmapping(df_subcell, df_mapping)
```

The first implementation may provide the join as a module-level helper rather
than hiding it inside the flat-file parser.

## Output Contract

Use one long table. A row represents:

```text
UniProt accession x subcellular location text x evidence
```

Columns:

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

Column semantics:

- `UniProtId`: accession for this row. Include secondary accessions so joins
  can match caller input when needed.
- `PrimaryUniProtId`: first accession from the `AC` block.
- `UniProtEntryName`: value from `ID`.
- `GeneName`: primary gene name from `GN   Name=...` when available.
- `ProteinName`: recommended full protein name from `DE   RecName: Full=...`
  when available.
- `SubcellularLocation`: location text from the subcellular location comment.
- `SubcellularLocationNote`: `Note=` text associated with the subcellular
  location comment when available.
- `EvidenceCode`: ECO code, for example `ECO:0000269`.
- `EvidenceSource`: evidence source namespace, for example `PubMed`.
- `EvidenceId`: evidence source identifier, for example a PMID.
- `SourceDb`: caller-supplied source label such as `sprot`.

If a location has multiple evidence references, emit multiple rows. If a
location has no evidence, emit one row with null evidence fields.

## Parsing Rules

Parse records delimited by `//`.

Entry identity fields:

- `ID`: entry name
- `AC`: all accessions; first accession is primary
- `GN`: primary gene name from `Name=`
- `DE   RecName: Full=`: recommended protein name

Subcellular location comments:

- collect continuation lines belonging to `CC   -!- SUBCELLULAR LOCATION:`
- keep parsing within the same comment topic until the next `CC   -!-` topic or
  record end
- split location statements conservatively on top-level periods only
- extract evidence blocks of the form `{ECO:...|Source:Id, ...}`
- extract `Note=` separately when present
- preserve unrecognized text as `SubcellularLocation` or
  `SubcellularLocationNote`; do not drop it silently

The parser should tolerate gzip and plain text through the existing
`open_uniprot_dat()` helper.

## Tidy Writer

Add schema constants:

```text
SCHEMA_VERSION_SUBCELLULAR_LOCATION = "uniprot-subcellular-location-v0.1"
COLS_SUBCELLULAR_LOCATION = [...]
SCHEMA_SUBCELLULAR_LOCATION = {...}
```

`write_subcellular_location_tidy()` should follow the existing UniProt flat-file
writer pattern:

- require a `.dat` snapshot
- write through a temporary TSV or directly build a LazyFrame if practical
- emit `data.parquet`
- support `should_write_manifest`
- support optional `should_hash_assets`
- include source `path`, media type, and source `bytes` when available

The manifest asset should be:

```json
{
  "path": "data.parquet",
  "kind": "canonical",
  "is_optional": false
}
```

## Relationship To GO And idmapping

This extractor is independent of GO xrefs:

- UniProt GO xrefs are `DR   GO; ...`
- subcellular location annotations are `CC   -!- SUBCELLULAR LOCATION: ...`

Do not make this a child of a GO xref API. If a later workflow needs a GO view,
add a separate mapping function that consumes this extracted table and an
explicit Location2GO resource.

The idmapping join is also a separate boundary. The flat-file parser should
produce accession-keyed annotations. A join helper can map those annotations to
protein IDs using an existing `idmapping_selected` extract.

## Goal-Mode Work Chain

This work can be split across parallel subagents after the design is accepted.

### Chain 0: Fixture And Real-Data Recon

Purpose: confirm real syntax before implementation hardens assumptions.

Inputs:

- local Swiss-Prot flat file path above
- existing tiny UniProt DAT fixture builder in
  `tests/integration/uniprot/test_uniprot.py`

Tasks:

- sample records with `SUBCELLULAR LOCATION`
- catalog common forms: simple location, multiple locations, evidence, `Note=`,
  isoform text, topology text
- produce 5-10 representative synthetic fixture records

Output:

- fixture examples for unit tests
- short notes on parser edge cases

Can run in parallel with Chain 1 after this document is available.

### Chain 1: Parser And Schema

Purpose: implement the core flat-file extraction.

Tasks:

- add constants for `COLS_SUBCELLULAR_LOCATION`,
  `SCHEMA_SUBCELLULAR_LOCATION`, and schema version
- add parser helpers for entry identity, subcellular location comments, and ECO
  evidence
- add `read_subcellular_location_frame()`
- add TSV writer and scanner only if needed for parquet writing
- keep existing eggNOG xref behavior unchanged

Output:

- tested parser functions in `src/bioextract/uniprot/util.py`
- focused unit tests for tiny fixtures

Depends on:

- Chain 0 examples for robust fixtures, but can start with minimal examples.

### Chain 2: Public API And Tidy Writer

Purpose: expose the extractor through `UniProtDatabase`.

Tasks:

- add `extract_subcellular_location()`
- add `write_subcellular_location_tidy()`
- add manifest wiring with `data.parquet`
- preserve `.dat` versus `idmapping` method guardrails
- optionally add `join_subcellular_location_idmapping()` as a module-level
  helper if the join contract is stable enough

Output:

- public API on `UniProtDatabase`
- tidy writer emitting `subcellular_location/data.parquet`

Depends on:

- Chain 1 parser and schema.

### Chain 3: Documentation And Examples

Purpose: make the feature usable without reading parser code.

Tasks:

- update `README.md`
- update `docs/architecture/uniprot-db.md`
- update `docs/testing/uniprot-db.md`
- document that this is Swiss-Prot curated annotation, not GO CC enrichment
- document output directory and `data.parquet`

Output:

- user-facing examples and test plan updates

Can run in parallel after Chain 2 API names are stable.

### Chain 4: Real-Data Smoke

Purpose: prove the extractor is useful on the real resource.

Tasks:

- run `from_dat(...).extract_subcellular_location()` on the real
  `uniprot_sprot.dat.gz`
- record row count, distinct primary accessions, evidence coverage, and a small
  deterministic preview
- run `write_subcellular_location_tidy()` into `/tmp`
- verify `data.parquet` schema and row count

Output:

- smoke-test evidence suitable for final handoff

Depends on:

- Chain 1 and Chain 2.

## Acceptance Criteria

Implementation is complete when:

- tiny fixture tests cover simple, multi-location, evidence, no-evidence, and
  note-bearing comments
- `.dat` and `.dat.gz` inputs both work
- `extract_subcellular_location()` returns the declared columns
- `write_subcellular_location_tidy()` writes `data.parquet` and optional
  manifest
- existing UniProt idmapping and eggNOG xref tests still pass
- real Swiss-Prot smoke succeeds on the 2026_01 local snapshot
- documentation states the scientific boundary and non-goals

## Risks

- UniProt comment grammar is richer than first fixtures. Mitigation: preserve
  source text conservatively and do not over-normalize.
- Swiss-Prot coverage is high quality but incomplete for some products.
  Mitigation: document coverage limits and keep missing annotations distinct
  from negative localization.
- Joining to project protein IDs may depend on caller-specific identifiers.
  Mitigation: keep accession-keyed output primary and make idmapping joins
  explicit.
