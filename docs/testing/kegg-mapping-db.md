# KEGG Mapping Test Standard

Version: v2.1
Date: 2026-08-14
Status: current

Changes to KEGG mapping access must verify:

1. both factory signatures, overlay precedence, direct-root semantics, and
   rejection of duplicate physical role files;
2. lazy construction and `collect_schema()` without organism enumeration or
   biological file reads;
3. bounded non-recursive directory discovery, strict rectangular organism
   profiles, selected-organism pruning, and empty-directory failure;
4. exact role column counts and identifier grammars with role/path/cause errors;
5. four-column gene-list parsing, aliases, opaque positions, duplicate
   summaries, and cross-organism rejection;
6. exact `organism`, `gene_annotation`, and `ko_annotation` schemas including
   every nested `List[Struct]` dtype and column position;
7. null/unavailable, empty/observed, and non-empty semantics for every
   capability combination;
8. direct gene pathways remain distinct from KO-mediated pathways;
9. namespace-aware normalization, all three selection namespaces, grouped
   lineage, multi-species matches, scoped unmatched IDs, and capability errors;
10. source/publication parity for all public LazyFrame relations and selections;
11. declarative publication source membership, nullable bytes/digest fields,
    no provenance-only content/stat pass, and no duplicate JSON inventory or
    biofetch dependency;
12. exact three-table and seven-capability publication inventories, row counts,
    nested schemas, organism scope mode, and read-only connections;
13. tamper rejection for metadata, capabilities, tables, schemas, counts,
    source rows, null invariants, and cross-species keys;
14. atomic `if_exists` behavior and staged-failure cleanup; and
15. rejection of the old wide `mapping` table/profile and absence of
    `mappings()`, mapping `build_tidy()`, limit, batch, and compatibility APIs.

Focused tests use temporary fixtures only. Real-scale validation is opt-in,
read-only on source data, writes only to a temporary destination, records wall
time, peak RSS, temporary/stage bytes, table counts, capabilities, and engine
settings, and verifies cleanup after both success and injected failure. Fixed
100- and 1,000-organism scopes are the optimization boundary; a complete source
build is a separately authorized release-readiness gate.
