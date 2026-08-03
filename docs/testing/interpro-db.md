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
- one capability-driven DuckDB publication
- DuckDB reopen, native read-only SQL, and domain-selection parity

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
- mapping-only construction writes only the `mapping` DuckDB relation.
- XML-capable construction writes `mapping`, `protein_term`, `term`, and
  `term_xref` in one DuckDB.
- `from_duckdb()` rejects non-v1 metadata and forged resource, profile,
  capability, source-role, table-role, schema/type, row-count, and column
  provenance inventories, including an incompatible persisted InterPro
  content-validation result.
- metadata v1 requires exact names, types, nullability, and primary keys for
  all five `_bioextract` tables.
- `connect()` returns distinct caller-owned read-only connections.
- a reopened handle rejects an atomically replaced publication path.
- relative publication paths remain bound across working-directory changes,
  and cached XML frames reject source identity changes.
- reopened selections preserve normalization, unique lookup, grouped fan-out,
  and unmatched behavior.
- `if_exists="fail"` and `"replace"` retain atomic publication behavior.
- retained lazy datasets reject changed source identities before atomic commit.
- invalid `namespace` raises targeted `ValueError`.
- duplicate positional Pfam matches collapse to one protein-term row.
- non-PFAM signatures are excluded.
- Pfam names remain distinct from InterPro entry names.
- missing or conflicting Pfam names, incomplete xrefs, malformed Pfam IDs, and
  mapping relationships absent from the XML raise targeted `ValueError`.
- a raw Pfam ID paired with the wrong InterPro ID fails even when both IDs exist
  independently in XML.
- requested source hashes cover every configured raw source.
- all three Pfam output frames are lazy and expose the exact public schemas.

## Real-Data Validation

The unified publication validates the explicitly assigned
`protein2ipr.dat.gz` and `interpro.xml.gz` roles by content. XML official
metadata supplies release identity; paths do not. Acceptance checks include
exact metadata, capability and schema inventories, readable row counts, a
deterministic sample after write, and:

- exact schemas for `mapping`, `protein_term`, `term`, and `term_xref`
- global `UniProtId + PfamId` uniqueness
- one non-empty name per published Pfam ID
- complete Pfam-to-InterPro xrefs
- equality with the PFAM projection of the full canonical mapping for a
  deterministic UniProt sample
- recorded output sizes, elapsed time, peak RSS, and selected-ID query time

Keep full-snapshot runs outside the default pytest suite. Historical canonical
and compact Parquet measurements, which predate unified DuckDB publication and
are not its size or reopen baseline, are recorded in the
[InterPro 108.0 benchmark](../benchmarks/20260714-v1.0-interpro-108-benchmark.md).
