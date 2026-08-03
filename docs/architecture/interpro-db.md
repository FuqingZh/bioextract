# InterProDatabase Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.interpro.InterProDatabase` provides path-first access to local InterPro
protein-domain mapping resources and exposes a stable UniProt-to-InterPro
annotation table for downstream enrichment and annotation joins.

The current version covers:

- `protein2ipr.dat.gz`
- optional `interpro.xml.gz` enrichment of entry type and member database
- full mapping extraction
- single and grouped UniProt selection
- capability-driven publication to one DuckDB file
- validated DuckDB reopening for domain selection and read-only SQL

It intentionally does not cover:

- InterProScan execution
- sequence-level domain search
- enrichment p-value calculation
- cross-ID conversion beyond UniProt accessions

## Raw Inputs

`InterProDatabase.from_mapping_files()` accepts exact file paths:

```text
protein2ipr.dat.gz
interpro.xml.gz
```

`protein2ipr.dat.gz` is required.

`interpro.xml.gz` is optional. When supplied, it contributes:

- `InterProType`
- `MemberDb`

If XML is absent, those columns remain null.

## Public API

```python
from bioextract import InterProDatabase

db = InterProDatabase.from_mapping_files(
    protein_to_interpro="protein2ipr.dat.gz",
    interpro_xml="interpro.xml.gz",
)

df_mapping = db.extract_mapping()

selection = db.select_ids(["P04637"], namespace="uniprot")
df_selected = selection.extract_mapping()
df_unmapped = selection.extract_unmatched_ids()
```

Grouped selections mirror the other DB contracts by prepending `GroupId`.

Publication writes every relation available from the constructed handle:

```python
from bioextract import InterProDatabase

db = InterProDatabase.from_mapping_files(
    protein_to_interpro="108.0/raw/protein2ipr.dat.gz",
    interpro_xml="108.0/raw/interpro.xml.gz",
)

result = db.write_duckdb("108.0/tidy/data.duckdb")

published = InterProDatabase.from_duckdb("108.0/tidy/data.duckdb")
with published.connect() as connection:
    table_names = connection.sql("SHOW TABLES").fetchall()
```

The arguments declare the two logical source roles. Paths and basenames do not
carry release identity. The official `INTERPRO` release declared in XML is
published as `release_version`, and mapping-to-XML relationships are validated
from content.

## Output Contract

The mapping table exposes:

```text
UniProtId
InterProId
InterProName
InterProType
MemberDb
MemberDbId
Start
End
```

The contract stays row-level:

- one UniProt accession can emit multiple InterPro entries
- one InterPro entry can emit multiple member-database rows
- positional coordinates remain attached to each row

## Publication

`write_duckdb(path)` publishes every relation supported by the constructed
source handle. A mapping-only handle publishes `mapping`. When InterPro XML is
also present, the same file additionally contains `protein_term`, `term`, and
`term_xref`. Missing XML produces absent Pfam relations, not misleading empty
tables and not a different container.

The write path is lazy:

1. scan `protein2ipr.dat.gz`
2. lazy-join XML-derived entry/member lookup frames when present
3. stream each relation through temporary transfer Parquet into staged DuckDB
4. validate metadata v1 and atomically commit the completed database

This keeps the main publication path aligned with the shared tidy contract
without forcing a full materialized DataFrame before write.

`build_tidy()` exposes the capability-driven lazy relation plan.

`from_duckdb()` accepts only metadata v1 publications with exact InterPro
resource identity, resource schema version, source profile, capability list,
source-role inventory, table names and roles, physical column types, row
counts, and column provenance. `connect()` returns a new caller-owned read-only
DuckDB connection on every call.

## Pfam Compact Contract

Pfam publication reads `protein2ipr.dat.gz` and `interpro.xml.gz` directly.
It does not consume or require a prior mapping publication.

`protein_term`:

```text
UniProtId
PfamId
```

`term`:

```text
PfamId
PfamName
```

`term_xref`:

```text
PfamId
InterProId
InterProName
InterProType
```

Only IDs matching `PF[0-9]{5}` and present in the raw protein mapping are
published. Positional repeats collapse to one `UniProtId + PfamId` row.
Pfam names come from the XML member `db_xref name`; InterPro names are kept
separately in the trace table. Missing names, conflicting names, incomplete
xrefs, invalid IDs, and cross-version inputs fail before asset writing. Raw
member relationships are validated against XML by the exact
`InterProId + PfamId` pair, so a valid Pfam ID cannot hide a mismatched
InterPro reference.

The publication resource schema is `interpro-v1`. The source profile and
`bioextract.capabilities` distinguish mapping-only from mapping-plus-Pfam
publications. Source hashes remain an optional writer cost.

All three published frames remain lazy. XML parsing materializes one compact
Pfam metadata index, and a streaming raw-file aggregation materializes one
compact used-pair index for validation and term filtering. Derived term and
xref outputs are not materialized before their Parquet sinks.

## Selection Contract

Accepted `namespace` values:

```text
uniprot
```

Single selection output prepends:

```text
InputId
InputNamespace
```

Grouped selection output prepends:

```text
GroupId
InputId
InputNamespace
```

`extract_unmatched_ids()` follows the same single/grouped shape as the
other DBs.

## Real Snapshot Status

The validated local InterPro snapshot is:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/interpro/mapping/108.0/raw
```

Published tidy output validation is recorded in the
[InterProDatabase test standard](../testing/interpro-db.md) and the
[InterPro 108.0 benchmark](../benchmarks/20260714-v1.0-interpro-108-benchmark.md).
