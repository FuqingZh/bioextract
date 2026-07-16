# InterProDb Architecture

Version: v1.0
Date: 2026-07-14
Status: current

## Goal

`bioextract.interpro.InterProDb` provides path-first access to local InterPro
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

`InterProDb.from_mapping_files()` accepts exact file paths:

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
from bioextract.interpro import InterProDb

db = InterProDb.from_mapping_files(
    file_protein2ipr="protein2ipr.dat.gz",
    file_interpro_xml="interpro.xml.gz",
)

df_mapping = db.extract_mapping()

selection = db.select_ids(["P04637"], kind_input_id="uniprot")
df_selected = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()
```

Grouped selections mirror the other DB contracts by prepending `GroupId`.

Compact Pfam publication is an InterPro tidy configuration:

```python
from bioextract.interpro import InterProDb

db = InterProDb.from_mapping_files(
    file_protein2ipr="108.0/raw/protein2ipr.dat.gz",
    file_interpro_xml="108.0/raw/interpro.xml.gz",
)

dataset = db.build_tidy(config="pfam")
report = db.write_tidy(
    "108.0/tidy/pfam",
    config="pfam",
    should_write_manifest=True,
    should_hash_sources=True,
    should_hash_assets=True,
)
```

Both inputs must be under the same `<version>/raw/` directory, and the
directory version must match the `INTERPRO` release declared in XML.

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

## Tidy Writing

`write_tidy()` emits:

```text
mapping.parquet
manifest.json
```

The write path is lazy:

1. scan `protein2ipr.dat.gz`
2. lazy-join XML-derived entry/member lookup frames when present
3. `sink_parquet()` directly to the final artifact

This keeps the main publication path aligned with the shared tidy contract
without forcing a full materialized DataFrame before write.

`build_tidy()` and `write_tidy()` default to `config="mapping"`.
`config="pfam"` selects the compact Pfam asset contract without adding a
parallel standalone API.

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

Accepted `kind_input_id` values:

```text
uniprot
```

Single selection output prepends:

```text
InputId
KindInputId
```

Grouped selection output prepends:

```text
GroupId
InputId
KindInputId
```

`extract_unmapped_input_ids()` follows the same single/grouped shape as the
other DBs.

## Real Snapshot Status

The validated local InterPro snapshot is:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/interpro/mapping/108.0/raw
```

Published tidy output validation is recorded in the
[InterProDb test standard](../testing/interpro-db.md) and the
[InterPro 108.0 benchmark](../benchmarks/20260714-v1.0-interpro-108-benchmark.md).
