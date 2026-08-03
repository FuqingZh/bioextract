# ChEBI Domain Access And Publication Plan

Storage status: superseded by the current
[Storage And Publication Convergence Plan](20260803-v1.0-storage-publication-convergence-implementation-plan.md).
The decisions and measurements below remain historical; metadata v1/v2/v3
statements and example filenames are not current publication authority.

Storage status: superseded by the current
[Storage And Publication Convergence Plan](20260803-v1.0-storage-publication-convergence-implementation-plan.md).
The decisions and measurements below remain historical; metadata v1/v2/v3
statements and example filenames are not current publication authority.

Date: 2026-07-29
Status: implemented
Authority:
[Domain Access Architecture](../architecture/20260729-v1.0-domain-access-architecture.md)

## Outcome

`ChEBIDatabase` will expose exact compound resolution and ChEBI-owned
relationships from a portable bioextract DuckDB. Callers will not need to join
primary IDs, secondary IDs, names, chemical properties, structures,
cross-references, or ontology edges themselves.

The first implementation targets the common annotation and ontology use cases:

- resolve primary and secondary ChEBI IDs to canonical compounds;
- retain caller and group lineage and report unmatched identifiers;
- extract stable compound facts, names, structures, and cross-references;
- extract direct ChEBI relations;
- traverse cycle-safe `is_a` ancestry and descendants;
- preserve star rating and obsolete-entry semantics.

Text search, chemical structure search, similarity calculation, enrichment, and
cross-database evidence integration remain outside this slice.

## Source Decision

### PostgreSQL dump is a build input, not a queryable publication

`pgsql_allstars.dump` is a PostgreSQL custom-format backup. It requires a
compatible `pg_restore` and a running PostgreSQL server before SQL queries can
be issued. It is not equivalent to SQLite or DuckDB and cannot be opened
read-only by `ChEBIDatabase`.

The inspected 2026-07-07 snapshot contains 15 normalized tables and useful
lookup relations, but it also contains internal or historical material that is
not automatically part of the public compound contract:

- 242,338 `compounds` rows versus 218,253 public flat-table compound rows;
- 14,354,566 `reference` rows;
- numeric `status_id`, `source_id`, and `relation_type_id` keys that must be
  joined to lookup tables;
- primary, secondary, merged, and parent-linked records that require explicit
  canonicalization rules.

Therefore the dump does not remove the need for a domain projection. Adding
`extract_*()` directly over raw restored tables would merely move repeated
joins and record-status decisions to every query.

### Canonical v1 build sources

The v1 domain publication uses:

1. ChEBI FULL OBO as the canonical self-describing compound and ontology
   source;
2. ChEBI FULL SDF only when molfile records are required;
3. official flat TSVs as compatible partial inputs and source-parity relations;
4. optional ChemOnt OBO as a separate ontology namespace.

The inspected FULL OBO contains 218,542 term stanzas, including 218,253
star-rated compound terms, 19,418 alternate IDs, 567,041 synonyms, 411,644
cross-references, 417,073 ontology edges, and chemical or structure property
values. It therefore avoids the unresolved numeric lookup identifiers present
in the standalone TSV files.

`from_release(source)` will discover the useful files that are actually
present. A release does not need every optional source. The resulting DuckDB
records capabilities and only enables extraction methods whose required
relations are present.

The PostgreSQL dump remains:

- a release audit and semantic verification source;
- a possible later source for proven user needs such as compound origins,
  comments, or the full reference relation;
- outside the v1 runtime and public constructor surface.

Do not add `from_postgresql_dump()` until a dump-exclusive, user-facing domain
need justifies the PostgreSQL restore dependency. If that need appears, design
the restore/import adapter as a separate optional build capability rather than
making PostgreSQL a dependency of publication-backed queries.

## Canonical DuckDB Relations

Biological relations live in `main`; publication provenance stays in
`_bioextract`.

| Relation | Role |
| --- | --- |
| `compound` | canonical ChEBI entity and scalar chemical facts |
| `secondary_id` | alternate or merged ChEBI ID to canonical ChEBI ID |
| `compound_name` | synonym or named representation with scope and source |
| `compound_cross_reference` | external source prefix and accession |
| `compound_relation` | directed, named ChEBI relation |
| `compound_structure` | optional SDF molfile records |
| `compound_wurcs` | optional WURCS representation |
| `chemont_term` | optional ChemOnt entity |
| `chemont_term_relation` | optional directed ChemOnt relation |
| `chemont_term_synonym` | optional ChemOnt synonym |
| `chemont_term_xref` | optional ChemOnt cross-reference |

ChEBI OBO is not duplicated into parallel `chebi_term*` and `compound*`
families. ChEBI is both the chemical database and the ontology; one canonical
compound identity and one named relation graph are less ambiguous.

OBO, SDF, and joined lookup fields are generated domain fields and use stable
`snake_case`. Official flat-table relations retained for partial source parity
keep their original headers, but public extraction never exposes unresolved
numeric source, status, or relation-type IDs as the only semantic value.

## Opening A Publication

```python
database = ChEBIDatabase.from_duckdb(path)
```

Opening verifies:

1. a readable DuckDB file and the required `_bioextract` metadata tables;
2. `resource_name = "chebi"` and supported metadata/schema versions;
3. actual `main` tables and views against the recorded inventory;
4. recorded row counts against actual row counts;
5. the set of domain capabilities derived from available relations.

Connections are short-lived and read-only. A publication-backed handle cannot
be passed back to `write_duckdb()`.

## Native DuckDB Access

Domain extraction is the recommended interface, not the exclusive interface.
`from_duckdb(path)` also exposes the complete publication through a native
read-only DuckDB connection:

```python
database = ChEBIDatabase.from_duckdb("chebi.duckdb")

with database.connect() as connection:
    frame = connection.sql(
        """
        SELECT compound.chebi_id, count(*) AS relation_count
        FROM compound
        LEFT JOIN compound_relation
          ON compound_relation.subject_chebi_id = compound.chebi_id
        GROUP BY compound.chebi_id
        """
    ).pl()
```

`connect()` is the advanced escape hatch and follows these rules:

- it is available only on a publication-backed handle;
- every call returns a new native `duckdb.DuckDBPyConnection`;
- the connection always opens the validated publication with
  `read_only=True`;
- the caller owns its lifetime and closes it directly or with `with`;
- all `main` relations, views, and `_bioextract` audit metadata remain
  queryable;
- no hidden shared connection is retained on `ChEBIDatabase`;
- concurrent callers create independent connections.

Do not add `ChEBIDatabase.sql()` or `ChEBIDatabase.query()`. Such methods would
only forward DuckDB calls while hiding connection and relation lifetimes.
Likewise, do not expose a persistent `.connection` attribute or a
`read_only=False` option.

A caller may deliberately open the file with `duckdb.connect(path)` for
write access, but after doing so bioextract no longer guarantees that domain
tables, `_bioextract.table_info`, and provenance remain consistent. The
validated `ChEBIDatabase` interface preserves publication immutability.

This creates three progressive levels:

1. `select_*()` and `extract_*()` for stable ChEBI semantics;
2. `connect()` for unrestricted read-only DuckDB SQL and relational queries;
3. direct `duckdb.connect(path)` for caller-owned mutation outside the
   bioextract publication contract.

## Selection Entry Points

```python
selection = database.select_compounds(
    ids,
    *,
    namespace,
    min_star_rating=1,
    include_obsolete=False,
)

selection = database.select_groups(
    ids_by_group,
    *,
    namespace,
    min_star_rating=1,
    include_obsolete=False,
)
```

`namespace` is explicit and one selection does not mix namespaces. Fixed
values are `chebi`, `inchi`, and `inchi_key`; external namespaces are the
normalized official prefixes present in `compound_cross_reference`.

Rules:

- `namespace="chebi"` accepts `CHEBI:15377` and bare numeric IDs, resolves
  primary and secondary IDs, and returns the canonical primary ID;
- `namespace="inchi"` and `"inchi_key"` perform exact matching;
- external cross-references use their prefix directly as `namespace`, for
  example `namespace="kegg.compound"` or `namespace="hmdb"`;
- an unknown prefix fails when the selection is created and reports all
  available values from that publication;
- names are not an identifier namespace because case, language, synonym scope,
  and non-uniqueness make name lookup a search problem;
- SMILES is not an identifier namespace because chemically equivalent strings
  require canonicalization rules.

`min_star_rating` is a ChEBI-owned curation filter and accepts `1`, `2`, or `3`.
The default retains the all-star publication. `include_obsolete=False` keeps
obsolete terms out of matches while secondary IDs still resolve to their
canonical active compound where the source declares that mapping.

## Deferred Selection And Eager Extraction

Selections are immutable query plans. They normalize and validate inputs but
do not query DuckDB until an `extract_*()` terminal is called. Terminals return
eager Polars `DataFrame` values for the same reasons documented by the Rhea
query contract: there is no genuine Polars database scan whose lifetime is
independent of the DuckDB connection.

Every matched output retains:

```text
GroupId          # grouped selections only
InputId
InputNamespace
ChEBIId
...
```

`extract_matches()` additionally returns how the input matched, including
primary ID, secondary ID, InChI, InChIKey, or cross-reference source.

## Extraction Surface

Core terminals:

```python
selection.extract_matches()
selection.extract_compounds()
selection.extract_names()
selection.extract_cross_references()
selection.extract_relations()
selection.extract_unmatched_ids()
```

`extract_compounds()` returns the one-row-per-selected-compound profile used
for annotation: canonical ID, preferred name, definition, star rating,
obsolete state, formula, charge, average mass, monoisotopic mass, SMILES,
InChI, InChIKey, and WURCS when available. This deliberately performs the
stable one-to-one joins callers would otherwise repeat.

Optional structure terminal:

```python
selection.extract_structures()
```

It returns all source structures and optional molfiles without inflating the
one-row compound profile.

Direct graph edges:

```python
selection.extract_relations(direction="both")
```

`direction` accepts `"outgoing"`, `"incoming"`, or `"both"`. The result always
retains explicit subject ID, relation type, and object ID; it does not relabel
all ChEBI predicates as parent-child relationships.

Cycle-safe ontology traversal:

```python
selection.extract_ancestors()
selection.extract_descendants()
```

Traversal is restricted to canonical `is_a` edges and tracks visited nodes so
malformed or cyclic source graphs cannot recurse forever. Other ChEBI
relations remain available through `extract_relations()`.

## Unmatched And Validation Contract

`extract_unmatched_ids()` returns a PascalCase `Reason` field with exactly:

```text
not_found
below_min_star_rating
obsolete_excluded
invalid_canonical_target
```

Any policy-valid candidate wins. Otherwise precedence is invalid target,
obsolete exclusion, star-rating exclusion, then no candidate.

Canonical compound IDs are fail-fast invariants. Missing or duplicate primary
IDs raise `bioextract.errors.IntegrityError`; staging files and WAL are removed
and an old destination is preserved. Orphan secondary IDs, cross-references,
relations,
structures, and WURCS rows are skipped and persisted as
`foreign_key_violation` warnings in `_bioextract.validation_issue`.

Metadata schema v3 always contains five tables and records
`bioextract.validation_status` plus
`bioextract.validation_issue_count`. Readers accept v1 as having no persisted
issue details, require `validation_issue` for v2, and reject unknown versions.

## Deliberate Non-goals

Do not add thin wrappers for ordinary projection, numeric filtering, sorting,
or arbitrary SQL.

Defer these features until they have a separate domain and dependency design:

- free-text or fuzzy name search;
- substructure, connectivity, or similarity search;
- SMILES canonicalization and chemical calculations;
- bulk export terminals for selections too large for memory;
- compound origins, comments, and the full dump reference history;
- ChEBI-to-Rhea or other cross-database application models.

The public ChEBI API exposes text, advanced, ontology, and structure search,
but faithfully reproducing those search semantics requires an indexed text
engine or chemistry toolkit. A misleading DataFrame substring or string-equality
wrapper is not an acceptable substitute.

## Implementation Sequence

1. Refactor the ChEBI OBO parser into canonical compound, alternate-ID, name,
   cross-reference, relation, property, and WURCS relations.
2. Add optional SDF parsing for molfile-backed `compound_structure`.
3. Make `from_release(source)` discover OBO, SDF, available TSVs, and optional
   ChemOnt without requiring a complete release.
4. Publish metadata v3 and integrity-check primary IDs,
   secondary-ID targets, relation endpoints, and structure ownership.
5. Add `ChEBIDatabase.from_duckdb()` with identity, schema, inventory, row
   count verification, and a fresh read-only native `connect()` escape hatch.
6. Implement exact compound and grouped selections with primary/secondary ID,
   InChI, InChIKey, and source-qualified cross-reference routing.
7. Implement core extraction terminals and capability errors.
8. Add cycle-safe direct and transitive relation queries.
9. Validate on the 2026-07-07 release and publish
   `tidy/chebi.duckdb`.
10. Rebuild Rhea v1 with complete `CHEBI:<number>` shared IDs and verify a
    direct equality join against `compound.chebi_id`.
11. Reassess dump-exclusive tables only after the portable v1 contract is
    measured against real downstream annotation use.

## Verification

Unit and integration tests must cover:

- primary, bare numeric, secondary, InChI, InChIKey, and source-qualified
  cross-reference matching;
- duplicate inputs, deterministic ordering, grouped isolation, and unmatched
  input accounting;
- one-row compound profiles and one-to-many names, structures, cross-references,
  and relations;
- `min_star_rating` and obsolete-entry behavior;
- relation direction and cycle-safe ancestor/descendant traversal;
- missing-capability errors for partial publications;
- wrong-resource, unsupported-schema, corrupt-inventory, and row-count
  rejection in `from_duckdb()`;
- `connect()` returning independent native read-only connections, supporting
  arbitrary SQL and Polars/Arrow result conversion, and rejecting persistent
  writes;
- source-backed handles rejecting `connect()` before publication;
- context-managed connections closing without leaving a shared handle;
- no query connection retained after extraction;
- failed publication preserving an existing destination and removing staging
  database and WAL files;
- real-snapshot checks for water (`CHEBI:15377`), one secondary ID, one external
  cross-reference, one `is_a` parent, and one descendant;
- `_bioextract` identity, source inventory, capability scope, and table counts;
- metadata v1/v2 compatibility, metadata v3 five-table integrity, validation
  warning persistence, and unknown-version rejection;
- a real-release read/extract smoke after publishing `tidy/chebi.duckdb`.
