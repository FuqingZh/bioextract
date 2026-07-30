# InterProDatabase Test Standard

Version: v1.0
Date: 2026-07-14
Status: current

## Scope

The InterProDatabase test standard covers:

- lightweight construction
- `protein2ipr.dat.gz` parsing
- optional XML enrichment
- full mapping extraction
- single and grouped UniProt selection
- unmapped reporting
- lazy tidy writing of one canonical parquet
- raw-to-Pfam compact generation without a full mapping prerequisite

It does not cover:

- InterProScan runtime behavior
- sequence search
- external service calls
- enrichment statistics

## Unit Tests

- `from_mapping_files()` accepts required and optional files.
- `extract_mapping()` returns the normalized mapping columns.
- XML enrichment fills `InterProType` and `MemberDb` when present.
- `select_ids(..., namespace="uniprot")` returns filtered rows.
- `select_groups(..., namespace="uniprot")` preserves `GroupId`.
- unmapped IDs are reported correctly.
- `write_parquet(path)` writes canonical mapping output with provenance.
- invalid `namespace` raises targeted `ValueError`.
- duplicate positional Pfam matches collapse to one protein-term row.
- non-PFAM signatures are excluded.
- Pfam names remain distinct from InterPro entry names.
- missing or conflicting Pfam names, incomplete xrefs, malformed Pfam IDs, and
  cross-version raw inputs raise targeted `ValueError`.
- a raw Pfam ID paired with the wrong InterPro ID fails even when both IDs exist
  independently in XML.
- formal manifests contain hashes for both raw sources and all three assets.
- `config="mapping"` remains the default and `config="pfam"` selects the
  compact contract.
- unknown configurations and Pfam requests without XML fail explicitly.
- all three Pfam output frames are lazy and expose the exact public schemas.

## Real-Data Validation

Canonical and compact publication are validated against same-version
`protein2ipr.dat.gz` and `interpro.xml.gz` files. Canonical acceptance checks
include exact schema, readable row count, and a deterministic sample after
write. Compact acceptance checks include:

- exact schemas for `protein_term.parquet`, `term.parquet`, and
  `term_xref.parquet`
- global `UniProtId + PfamId` uniqueness
- one non-empty name per published Pfam ID
- complete Pfam-to-InterPro xrefs
- equality with the PFAM projection of the full canonical mapping for a
  deterministic UniProt sample
- recorded output sizes, elapsed time, peak RSS, and selected-ID query time

Keep full-snapshot runs outside the default pytest suite. The observed canonical
and compact 108.0 baselines are recorded in the
[InterPro 108.0 benchmark](../benchmarks/20260714-v1.0-interpro-108-benchmark.md).
