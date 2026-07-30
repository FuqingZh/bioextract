# RheaDatabase Architecture

Version: v1
Date: 2026-07-29
Status: current

## Purpose

`RheaDatabase` exposes official local Rhea releases as one self-contained,
query-ready DuckDB database. It encodes reaction direction, participants,
compounds, publications, and resource-owned cross-references so callers do not
repeat RDF/XML, SDF, TSV, gzip, and internal-join logic.

DuckDB is the publication container rather than the parsing engine. RDF/XML and
SDF are streamed with Python readers; large tabular mappings are loaded directly
by DuckDB, including gzip-compressed TrEMBL input. DuckDB retains columnar
compression and vectorized query execution while avoiding a directory of
loosely related Parquet files.

Normalized RDF/SDF batches use short-lived local Parquet transfer files so
DuckDB can ingest them column-wise without a PyArrow runtime dependency or slow
row-by-row SQL. These files are deleted immediately and are not publication
artifacts; the only output is the DuckDB database.

## Progressive Entry Points

The public constructors disclose complexity by capability:

- `from_reaction_files(...)` requires RDF and the direction quartet; hierarchy,
  obsolete IDs, and reaction SMILES are optional.
- `from_compound_files(...)` accepts any non-empty combination of SDF,
  ChEBI names, and the pH 7.3 mapping.
- `from_cross_reference_files(...)` accepts any non-empty combination of the
  aggregate xref, Swiss-Prot, and TrEMBL mappings.
- `from_release(source)` accepts an extracted directory or zip/tar archive and
  strictly requires all 15 assets in the current official release contract.
- `from_duckdb(path)` validates a bioextract Rhea publication and opens it as a
  read-only domain-query handle.

Partial inputs create only meaningful tables. A missing component is represented
by an absent table, not an empty table that could be mistaken for an observed
zero-row dataset.

Compression is detected internally. Partial constructors therefore do not
encode `.gz` in their method names, and release archives need not be manually
extracted.

## Naming Contract

Python uses conventional snake_case method and parameter names. DuckDB schemas,
tables, columns, and views also use snake_case. Original source headers remain
source provenance rather than becoming the public query contract; this avoids
leaking inconsistent uppercase and resource-specific conventions into user SQL.

Stable biological identifiers keep their domain names (`rhea_id`, `chebi_id`,
`go_id`, `uniprot_id`). No implicit identifier conversion is performed.

## Database Inventory

Tables are conditional on their source:

| Source | Tables |
| --- | --- |
| Rhea RDF/XML | `reaction`, `reaction_side`, `compound`, `compound_reactive_part`, `reaction_participant`, `reaction_publication` |
| Direction TSV | `reaction_quartet` |
| Relationship TSV | `reaction_relationship` |
| Obsolete TSV | `obsolete_reaction` |
| Reaction SMILES TSV | `reaction_smiles` |
| Rhea SDF | `compound_structure` |
| ChEBI name TSV | `chebi_name` |
| ChEBI pH 7.3 TSV | `chebi_ph7_3_mapping` |
| Aggregate xref TSV | `reaction_xref` plus `reaction_ec` and `reaction_go` views |
| UniProt TSVs | `reaction_uniprot` |

When RDF is present, `reaction_participant_direction` retains exact reaction,
master, direction, side, participant, and compound-link fields. Its nullable
`directional_role` avoids forcing callers to repeatedly join exact reaction
direction to master-reaction sides.

List-valued source facts are normalized only when they represent independent
many-to-many entities, such as publications and participants. Scalar chemical
fields stay on their owning record. Symbolic coefficients and charges retain
their text forms; numeric companion columns are populated only when conversion
is lossless.

## Direction And Hierarchy

Rhea is not one tree. Reaction relationships form a direction-aware directed
acyclic graph, while each master reaction also has an exact four-member
direction family (`UN`, `LR`, `RL`, `BI`).

RDF side references are authoritative for deriving each active exact reaction's
`master_id` and `direction`. Historical obsolete tombstones without side
references retain null values rather than receiving an invented direction.
When the direction TSV is also supplied, every current quartet is checked
against the RDF-derived semantics. Hierarchy edges remain explicit in
`reaction_relationship`; no single-parent tree projection is invented.

Side labels alone do not define substrate and product. The
`reaction_participant_direction` view applies:

| Exact direction | Left side | Right side |
| --- | --- | --- |
| `LR` | `substrate` | `product` |
| `RL` | `product` | `substrate` |
| `UN` | null | null |
| `BI` | null | null |

Exact Rhea IDs are never collapsed to their master reaction during selection
or participant extraction.

## Publication Query Interface

`RheaDatabase.from_duckdb(path)` accepts only a supported bioextract Rhea
publication. Opening validates the metadata-version-specific `_bioextract`
tables, resource and
schema identity, the complete `main` table inventory, and every recorded row
count. Partial publications remain valid; operations fail with
`RheaCapabilityError` only when their required relations are absent.

The two selection entry points are:

```python
database.select_reactions(
    ids,
    namespace="chebi",
    include_obsolete=False,
)
database.select_groups(
    ids_by_group,
    namespace="ec",
    include_obsolete=False,
)
```

One selection uses one explicit namespace. Supported namespaces are `rhea`,
`chebi`, `uniprot`, `ec`, `go`, `ecocyc`, `kegg_reaction`, `macie`,
`metacyc`, and `reactome`. ChEBI lookup is exact against participant compound
identity; the pH 7.3 mapping is not applied implicitly.

Selections are immutable deferred query plans. They normalize identifiers,
validate capabilities, and retain query policy without executing DuckDB.
Terminal methods execute through short-lived read-only connections and return
eager Polars `DataFrame` objects:

- `extract_matches()`;
- `extract_reactions()`;
- `extract_participants()`;
- `extract_cross_references()`;
- `extract_unmatched_ids()`;
- `extract_publications()`;
- `extract_relationships()`.

All matched outputs retain `InputId`, `InputNamespace`, `RheaId`, `MasterId`,
and `Direction`; grouped selections additionally retain `GroupId`. There is no
`scan_*()` facade because Polars does not provide a true database scan whose
plan can safely outlive the short-lived DuckDB connection. A future
bounded-memory contract should use a real iterator or direct writer terminal.

## Complete Release Validation

The complete release contract currently consists of:

`LICENSE.txt`, `chebiId_name.tsv`, `chebi_pH7_3_mapping.tsv`,
`rhea-directions.tsv`, `rhea-obsoletes.tsv`, `rhea-reaction-smiles.tsv`,
`rhea-relationships.tsv`, `rhea-release.properties`, `rhea.rdf`, `rhea.sdf`,
`rhea2ec.tsv`, `rhea2go.tsv`, `rhea2uniprot_sprot.tsv`,
`rhea2uniprot_trembl.tsv`, and `rhea2xrefs.tsv`.

Plain and gzip variants are accepted where applicable. Duplicate logical assets
and incomplete releases are rejected. `rhea2ec` and `rhea2go` validate the
matching rows in `rhea2xrefs`; they do not create duplicate physical tables.

## Publication And Provenance

`write_duckdb(...)` builds and validates a separate staging database, closes
DuckDB, and then atomically replaces the destination. Table-sized transactions
allow large mappings to flush without retaining an entire release in one
transaction. `if_exists="fail"` is the default; `"replace"` never destroys the
previous database before a successful replacement exists.

The `_bioextract` schema contains:

- `metadata`: metadata schema version, resource name, resource schema version,
  scope, generation time, package version, and release metadata when
  available;
- `source_file`: logical source name, display path, byte size, media type, and
  optional SHA-256;
- `table_info`: physical table name, semantic role, and row count;
- `column_mapping`: necessary source-to-output column mappings; the table is
  present and empty when no mapping was required.
- `validation_issue`: non-fatal source-integrity findings; Rhea v1 creates it
  empty because current integrity failures remain fail-fast.

New publications use metadata schema v2. Readers accept metadata v1 without
persisted issue details and v2 with all five tables; unknown versions fail.

Rhea v1 stores `compound.chebi_id` and `underlying_chebi_id` as complete
`CHEBI:<number>` CURIE strings. The resource schema version remains
`rhea-duckdb-v1`; an older numeric physical layout is rejected because it no
longer satisfies the v1 shared-ID contract.

`_bioextract` is an application-owned internal provenance schema. It is not a
DuckDB or SQL standard namespace and contains no Rhea biological relations.

Source hashing is opt-in because the TrEMBL mapping is large and hashing adds a
full extra read.

Metadata keys use the shared `bioextract.*` namespace, including
`bioextract.metadata_schema_version`, `bioextract.resource_name`,
`bioextract.schema_version`, `bioextract.scope`, and source/release identity.
A publication-backed handle cannot be republished; rebuilding must begin from
official source files.

## Non-goals

`RheaDatabase` does not download releases, call Rhea services, calculate
reaction networks, infer missing chemistry, merge ChEBI ontology content, or
impose an application-specific schema. ChEBI/ChemOnt publication is a separate
resource boundary even when a caller later joins it by `chebi_id`.
