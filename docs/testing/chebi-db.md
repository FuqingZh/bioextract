# ChEBIDatabase Test Standard

Version: v1
Date: 2026-07-29
Status: current

Automated tests cover:

- representation-specific OBO and table constructors, directory/archive
  discovery, explicit role overlays, candidate ambiguity, provenance, plain
  and gzip inputs, SDF molfile supplementation, and unsafe-member rejection;
- canonical compound, secondary ID, name, xref, relation, structure, WURCS,
  and separate ChemOnt relations;
- complete `CHEBI:<number>` keys and dynamic external namespaces;
- primary/secondary, InChI, InChIKey, grouped selection, star rating,
  obsolete policy, and all four unmatched `reason` values;
- direct relations plus cycle-safe `is_a` ancestor/descendant traversal;
- canonical fail-fast behavior, old-target preservation, staging/WAL cleanup,
  orphan-row skipping, and validation issue persistence;
- metadata v2 exact five-table requirements, rejection of every other version,
  inventory validation, and count validation;
- independent native read-only `connect()` calls and arbitrary SQL;
- source-backed handles rejecting query operations before publication.

Before a release is accepted, build the official FULL OBO plus SDF into
the versioned `tidy/data.duckdb`, open it with `from_duckdb()`, query water
(`CHEBI:15377`), exercise a secondary ID, external prefix, relation traversal,
and unmatched case, then compare every `_bioextract.table_info` count with its
live table. Finally verify an equality-only Rhea participant-to-ChEBI compound
join with no casts or prefix construction.
