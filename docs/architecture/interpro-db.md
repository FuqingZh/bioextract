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
- flat tidy writing to one canonical `mapping.parquet`
- direct raw-to-Pfam compact publication

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

Compact Pfam publication is an InterPro tidy configuration:

```python
from bioextract import InterProDatabase

db = InterProDatabase.from_mapping_files(
    protein_to_interpro="108.0/raw/protein2ipr.dat.gz",
    interpro_xml="108.0/raw/interpro.xml.gz",
)

result = db.write_duckdb(
    "108.0/interpro_pfam.duckdb",
    config="pfam",
)
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

`write_parquet(path)` publishes the independent InterPro mapping.
`write_duckdb(path, config="pfam")` publishes the related Pfam
`protein_term`, `term`, and `term_xref` relations together.

The write path is lazy:

1. scan `protein2ipr.dat.gz`
2. lazy-join XML-derived entry/member lookup frames when present
3. `sink_parquet()` to a staging artifact before atomic publication

This keeps the main publication path aligned with the shared tidy contract
without forcing a full materialized DataFrame before write.

`build_tidy()` exposes the lazy relation plan; `write_parquet(path)` and
`write_duckdb(path)` publish its single- and multi-relation products.

## Pfam Compact Contract

Pfam publication reads `protein2ipr.dat.gz` and `interpro.xml.gz` directly.
It does not consume or require the full `tidy/mapping.parquet`.

`protein_term.parquet`:

```text
UniProtId
PfamId
```

`term.parquet`:

```text
PfamId
PfamName
```

`term_xref.parquet`:

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

The compact manifest schema is `interpro-pfam-v0.1`. Source and asset hashes
are optional API costs and are enabled for formal resource publication.

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
