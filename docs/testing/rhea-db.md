# RheaDatabase Test Standard

Version: v1
Status: current

Changes to `RheaDatabase` must verify:

1. `from_files()` has exactly 11 optional keyword-only roles, rejects an empty
   profile and duplicate physical files across roles, requires RDF and
   directions for any reaction role, and accepts independent compound and
   cross-reference roles plus mixed partial profiles.
2. gzip/plain RDF detection is content-based.
3. RDF direction, participants, coefficients, locations, compounds, symbolic
   charges, and direction-aware roles survive normalization; `UN` and `BI`
   roles are null.
4. A partial write creates no tables for absent components.
5. Aggregate xrefs create EC/GO views and UniProt provenance distinguishes
   Swiss-Prot from TrEMBL.
6. A complete extracted release and its archive produce the same row counts.
7. Complete-release discovery rejects missing or duplicate logical assets.
8. Existing output follows `fail`/`replace`, and a failed staging build does not
   corrupt the destination.
9. Source hashes are absent by default and valid when requested.
10. `from_duckdb()` rejects wrong resource identity, incompatible schema,
    missing metadata tables, inventory drift, and row-count drift.
11. Partial publications open successfully, while unsupported namespace or
    extraction operations raise `RheaCapabilityError`.
12. Rhea, ChEBI, UniProt, and external-xref selections preserve exact reaction
    ID, master ID, direction, input namespace, unmatched IDs, and group
    isolation.
13. ChEBI selection performs no implicit pH 7.3 conversion.
14. Domain extraction agrees row-for-row with direct DuckDB SQL for reaction,
    participant, cross-reference, publication, and hierarchy relations.
15. All public classes and methods provide direct, observable examples.
16. Rhea ChEBI fields are complete CURIE strings and numeric legacy
    publications are rejected as an incompatible v1 physical layout.
17. Metadata v1 remains readable, metadata v2/v3 requires the fifth
    `validation_issue` table, and unknown metadata versions are rejected.
18. `connect()` returns independent native read-only connections and rejects
    persistent writes.
19. The three superseded explicit-file constructors are absent, all supplied
    roles appear in provenance with `obsolete_reactions` recorded as
    `obsoletes`, and each source profile has `reactions`, `compounds`,
    `cross_references`, or `partial` construction scope. `write_duckdb()`
    embeds that value as `bioextract.scope`; `from_duckdb()` returns a handle
    whose snapshot scope is `publication` while retaining the embedded
    construction scope for audit.

Run the focused suite with:

```bash
PYTHONPATH=src pytest tests/test_rhea.py tests/test_docstring_examples.py
```

Before publication, also run the full project checks and a real-snapshot smoke
test. The smoke test should build a database from the current complete release,
open it through `RheaDatabase.from_duckdb()`, compare
`_bioextract.table_info` to live counts, and exercise at least one Rhea,
ChEBI, UniProt, EC, participant, and unmatched-ID query.
The release smoke must also directly join participant `chebi_id` values to a
ChEBI `compound.chebi_id` column without casts or prefix concatenation.
