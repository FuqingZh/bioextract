# KEGG Metabolic Domain Access Plan

Storage naming status: superseded by the current
[Storage And Publication Convergence Plan](20260803-v1.0-storage-publication-convergence-implementation-plan.md).
The decisions, artifact paths, and measurements below remain historical;
resource-named DuckDB guidance is not current schema or filename authority.

Storage naming status: superseded by the current
[Storage And Publication Convergence Plan](20260803-v1.0-storage-publication-convergence-implementation-plan.md).
The decisions, artifact paths, and measurements below remain historical;
resource-named DuckDB guidance is not current schema or filename authority.

Date: 2026-07-30
Status: completed
Authority:
[Domain Access Architecture](../architecture/20260729-v1.0-domain-access-architecture.md)
and [Materialized Dataset Contract](../architecture/tidy-dataset-contract.md)

## Outcome

`KEGGDatabase` will turn one local KEGG metabolic snapshot into a portable,
queryable domain publication. Callers should be able to start from a KEGG,
ChEBI, PubChem, Rhea, EC, KO, module, or reference-pathway identifier and
retrieve the KEGG-owned relationships without parsing flat records or
repeating compound-reaction-enzyme-module joins.

The canonical full publication is:

```text
kegg.duckdb
```

It remains separate from KEGG BRITE and organism mapping products. Those
products have independent release scopes and flat analytical contracts:

```text
kegg_brite.parquet
kegg_gene_annotation.parquet
kegg.duckdb
```

Here `kegg.duckdb` means the publication inside the `kegg/metabolic/<version>`
product directory. It is not a universal file that absorbs independently
versioned BRITE, organism mapping, pathway, or other KEGG products. A caller
may attach or join those publications explicitly when it needs a cross-product
query.

### Filename convention

The managed resource layout already carries product and release identity:

```text
kegg/metabolic/2026-07/tidy/kegg.duckdb
```

Within this layout, durable DuckDB publications use `<resource>.duckdb`.
The product qualifier is not repeated in the filename because `metabolic` is
already an unambiguous directory component. `data.duckdb` may be used for an
opaque temporary artifact, but not as the recommended durable publication
name: a copied or downloaded `<resource>.duckdb` remains recognizable without
its original parent directories.

Add the smallest product qualifier only if multiple canonical DuckDB files
must coexist in the same `tidy` directory. This is a deployment convention,
not a writer default or compatibility identifier; `write_duckdb(path)` still
requires the caller to provide the complete destination.

## Inspected 2026-07 Snapshot

The design is based on:

```text
/cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/kegg/metabolic/2026-07
```

The snapshot contains:

| Entry family | Entries | Batch files |
| --- | ---: | ---: |
| compound | 19,619 | 1,962 |
| reaction | 12,459 | 1,246 |
| enzyme | 8,343 | 835 |
| module | 573 | 58 |

It also contains seven global relation files:

| Source relation | Rows |
| --- | ---: |
| `compound_pubchem.tsv` | 19,443 |
| `compound_reaction.tsv` | 52,626 |
| `reaction_enzyme.tsv` | 10,752 |
| `reaction_ko.tsv` | 12,276 |
| `reaction_module.tsv` | 3,187 |
| `reaction_pathway.tsv` | 19,828 |
| `module_pathway.tsv` | 1,681 |

Each entry batch contains at most ten KEGG flat-file records. `manifest.lock`
records the release identity, source URLs, byte sizes, and SHA-256 values, but
it is optional provenance input rather than a required biofetch contract.

## Domain Decisions

### KEGG reactions do not provide physiological direction

All 12,459 inspected `EQUATION` fields use `<=>`. This agrees with the official
KEGG REACTION definition: a KEGG reaction is assumed reversible and its two
sides are separated by `<=>`.

The publication therefore records `left` and `right` sides. It must not call
them substrate/product roles or infer physiological direction. Directional
module steps and directional Rhea accessions are separate concepts.

See the official
[KEGG REACTION entry help](https://www.kegg.jp/kegg/document/help_bget_reaction.html).

### Module definitions are logical expressions

KEGG module `DEFINITION` is not a flat KO list. The official syntax uses:

- a top-level space for connected required blocks;
- `+` for components of a molecular complex;
- `,` for alternatives;
- `-` for an optional component;
- parentheses for grouping;
- K or, in some signature modules, M identifiers as leaves.

KEGG treats space and plus as AND-like operations and comma as OR, while
preserving their different biological roles. A completeness check operates on
top-level blocks. See the official
[KEGG MODULE description](https://www.kegg.jp/kegg/module.html).

The publication keeps the original expression and a lossless parsed tree.
A flattened `module_ko` relation is useful for discovery but is not sufficient
for completeness evaluation.

### Global relations use official link files

The seven TSV link/conv files are the canonical source for the corresponding
many-to-many relations. Duplicate fields embedded in entry records are used
for validation or attributes; they are not published as a second competing
relationship.

Relations not available as global TSVs may be derived from the relevant entry
field when their semantics are explicit, including module definition members,
module reaction steps, compound/reaction external references, and enzyme-KO
links.

## Public Construction

Keep `KEGGDatabase` as the single top-level KEGG resource type.

```python
database = KEGGDatabase.from_metabolic_release(
    source,
    release_version="2026-07",  # optional caller-known official identity
)
```

`from_metabolic_release(source)` accepts a release directory, its `raw`
directory, or a supported archive containing that layout. It discovers
available logical roles rather than depending on an organization-specific
absolute path. `source` declares where to discover the logical roles; its
directory, file, or archive name never supplies release identity. When known,
the caller may pass `release_version`; otherwise publication omits release
metadata.

An explicit partial-input constructor supports tests, nonstandard layouts, and
callers that possess only some official assets:

```python
database = KEGGDatabase.from_metabolic_files(
    compound_list=...,
    compound_entries=...,
    reaction_list=...,
    reaction_entries=...,
    enzyme_list=...,
    enzyme_entries=...,
    module_list=...,
    module_entries=...,
    compound_pubchem=...,
    compound_reaction=...,
    reaction_enzyme=...,
    reaction_ko=...,
    reaction_module=...,
    reaction_pathway=...,
    module_pathway=...,
    release_version=...,
)
```

Each `*_entries` value describes one logical batch collection and may resolve
to a directory, archive, single batch, or sequence of batches. List and
relation parameters are explicit official file roles. `release_version` lets a
nonstandard partial layout retain caller-known source identity. Do not add
compression-specific constructors or four parallel
`from_compound_*`/`from_reaction_*` resource types.

A partial input still publishes a DuckDB. Metadata records its capabilities,
and domain methods raise a targeted capability error when a required relation
was not supplied.

## Opening A Publication

```python
database = KEGGDatabase.from_duckdb("kegg.duckdb")
```

Opening verifies resource identity, supported metadata and resource schema
versions, actual relation inventory, and row-count parity. It rejects a KEGG
BRITE or mapping Parquet and a DuckDB belonging to another resource.

As with ChEBI and Rhea:

```python
with database.connect() as connection:
    relation = connection.sql(
        """
        SELECT reaction.reaction_id, count(*) AS participant_count
        FROM reaction
        JOIN reaction_participant USING (reaction_id)
        GROUP BY reaction.reaction_id
        """
    )
```

`connect()` is publication-only, creates a new native DuckDB connection with
`read_only=True`, and leaves connection lifetime to the caller. Do not add
forwarding `sql()`/`query()` methods, a shared connection, or a writable flag.

This provides progressive disclosure:

1. `select_ids()` and domain `extract_*()` methods;
2. native read-only DuckDB relations through `connect()`;
3. direct caller-owned DuckDB access outside bioextract guarantees.

## Canonical Relations

Biological relations live in `main`; provenance remains in `_bioextract`.

### Entities and attached facts

| Relation | Role |
| --- | --- |
| `compound` | one KEGG C-number with primary name, formula, exact mass, and molecular weight |
| `compound_name` | ordered accepted and alternative names |
| `compound_cross_reference` | external or KEGG-internal reference with relationship semantics |
| `reaction` | one R-number with name, definition, raw equation, and reversible-as-defined flag |
| `reaction_name` | ordered reaction names |
| `reaction_cross_reference` | Rhea and other external reaction references |
| `reaction_class` | reaction-to-RCLASS relation and optional compound-pair text |
| `enzyme` | one EC number with status, class, systematic name, comment, and history |
| `enzyme_name` | accepted and alternative enzyme names |
| `enzyme_cross_reference` | external enzyme nomenclature references |
| `enzyme_replacement` | transferred or obsolete EC number to replacement EC number |
| `module` | one M-number with type, name, class, raw definition, and diagram identifier |
| `module_definition_node` | lossless ordered expression tree for the module definition |

Do not publish KEGG `ATOM`/`BOND`, ENZYME `GENES`/`REFERENCE`, MODULE
`COMPLETE`, or embedded BRITE trees in v1:

- KCF structure graphs need a separate demonstrated structure-search use case;
- ENZYME `GENES` and MODULE `COMPLETE` dominate input size and overlap the
  organism mapping product;
- bibliography is not needed for the first annotation contract;
- BRITE already has a separate product.

The streaming parser skips these fields without buffering their contents.

### Domain relations

| Relation | Role |
| --- | --- |
| `reaction_participant` | ordered equation participant, side, namespace, and symbolic or numeric coefficient |
| `reaction_enzyme` | R-number to EC number |
| `reaction_ko` | R-number to K-number |
| `reaction_module` | R-number to M-number membership |
| `reaction_pathway` | R-number to reference map ID |
| `module_pathway` | M-number to reference map ID |
| `module_ko` | flattened module member index for discovery, not evaluation |
| `module_reaction_step` | ordered, module-specific directional reaction step |
| `module_compound` | compound membership explicitly listed by a module |
| `enzyme_ko` | EC number to K-number from ENZYME orthology |

Do not create both `reaction_module` and `module_reaction` physical tables.
Direction of lookup is a query concern, not a second biological relation.
Likewise, `reaction_compound` is a view or domain projection over
`reaction_participant`, not a duplicate base table.

Pathway, KO, glycan, drug, and RCLASS identifiers may appear as endpoints
without a local entity table because their complete entry families are not in
this snapshot. A foreign identifier is not promoted into a stub entity that
pretends to contain metadata.

## Stable Identifiers

Generated columns use singular `snake_case`. IDs stay strings:

```text
compound_id   C00002
reaction_id   R00002
ec_number     1.1.1.1
module_id     M00001
ko_id         K00001
pathway_id    map00010
glycan_id     G00092
```

Never parse an EC number or symbolic stoichiometric coefficient as a floating
point value.

Cross-database shared IDs are normalized at publication time:

```text
ChEBI 15422  -> CHEBI:15422
RHEA 22455   -> RHEA:22455
```

This permits direct equality joins to current ChEBI `compound.chebi_id` and
Rhea `reaction.accession`. Namespace and identifier remain separate columns
in cross-reference relations; normalization does not create a universal
cross-database entity model.

## Reaction Participant Contract

`reaction_participant` has:

```text
reaction_id
side
position
participant_namespace
participant_id
coefficient_text
coefficient_numeric
```

`side` is `left` or `right`. `participant_namespace` is at least
`kegg_compound` or `kegg_glycan`; both C and G identifiers occur in the
inspected equations.

`coefficient_text` preserves values such as `n`, `(n+1)`, and `(n-1)`.
`coefficient_numeric` is populated only when the source coefficient is an
ordinary number. A missing coefficient means one. The parser never evaluates
symbolic expressions.

The raw equation remains on `reaction` so parsing can be audited. The official
`compound_reaction.tsv` set is compared with parsed C-number participants.
Differences become validation issues unless they demonstrate that the
participant contract itself is invalid.

## Module Definition Contract

`module_definition_node` represents one ordered abstract syntax tree:

```text
module_id
node_id
parent_node_id
position
node_kind
member_namespace
member_id
```

`node_kind` is:

```text
sequence
complex
alternative
optional
identifier
```

Leaves use `member_namespace=ko` or `module`. Parentheses affect the tree but
do not create a fake biological entity.

`evaluate_modules(ko_ids)` is allowed because it evaluates KEGG-owned module
logic rather than performing statistical enrichment. It returns exact
top-level-block evaluation:

```text
ModuleId
RequiredBlockCount
SatisfiedBlockCount
IsComplete
MissingBlockIndexes
```

Do not add an arbitrary completeness threshold or a custom stability score.
Callers can filter the returned exact facts themselves.

## Selection And Domain Extraction

`select_ids()` and `select_groups()` remain the common KEGG selection verbs.
On a metabolic publication they accept:

```text
kegg_compound
chebi
pubchem
kegg_reaction
rhea
ec
ko
kegg_module
kegg_pathway
```

Unknown namespaces report the values actually present in the publication.

The resulting deferred selection exposes:

```text
extract_matches()
extract_compounds()
extract_reactions()
extract_participants()
extract_enzymes()
extract_kos()
extract_modules()
extract_pathway_memberships()
extract_cross_references()
extract_unmatched_ids()
```

Selection is reaction-centered but not an unbounded graph closure:

1. input IDs resolve to canonical anchor entities;
2. reactions directly linked to those anchors form the selected reaction set;
3. participants, enzymes, KOs, modules, and pathways are extracted only from
   that selected reaction set;
4. compound anchors remain identifiable separately from co-participants;
5. every output retains `GroupId`, `InputId`, and `InputNamespace` lineage when
   applicable.

This rule supports the real annotation paths without a hidden recursive walk:

```text
ChEBI/PubChem/C number -> compound -> reactions -> EC/KO/module/pathway
Rhea/R number         -> reaction -> participants and functional links
EC/KO/M number        -> reactions -> compounds and pathways
map ID                -> reactions -> compounds, EC, KO, and modules
```

Primary KEGG identifiers match exactly after prefix normalization. External
cross-references may map one input to multiple canonical KEGG entities.
Transferred EC entries resolve to their valid replacement entries and record
`MatchType=replacement`. Deleted or obsolete entries without an accepted
target are reported through the existing bioextract reason vocabulary:

```text
not_found
obsolete_excluded
invalid_canonical_target
```

`include_obsolete=False` is the default selection policy. Advanced callers
that need the exact historical entry can opt in or query it directly through
DuckDB.

Selections are deferred query plans; `extract_*()` methods are eager terminal
operations returning Polars `DataFrame`. Large or custom lazy queries use the
native DuckDB relation returned by `connect().sql()`.

## Parsing And Build Strategy

1. Resolve input roles and sort batch paths deterministically.
2. Stream KEGG flat files using their fixed-width field labels and `///`
   record terminator.
3. Buffer only fields used by the canonical contract; skip large excluded
   fields while scanning.
4. Read global link/conv TSVs through lazy Polars scans.
5. Normalize identifiers and build canonical relations.
6. Validate entity uniqueness, list-to-entry parity, relation endpoints,
   equation participants, module grammar, and duplicated link fields.
7. Write a staging DuckDB adjacent to the destination.
8. Create `_bioextract` metadata schema v3, close and reopen read-only, verify
   inventory and counts, then atomically replace the target.

Tables are sorted by their primary lookup keys so DuckDB zone maps can prune
common selections. Do not add indexes until a real-snapshot benchmark shows a
point-query benefit that justifies their size.

## Integrity Policy

Fail the full canonical build for:

- missing or duplicate primary IDs;
- complete-release list/entry mismatch;
- an unparseable reaction equation needed for participant semantics;
- an unparseable module definition needed for exact evaluation;
- unsupported metadata or resource schema versions;
- table inventory or recorded row-count mismatch after writing.

For a dependent link whose endpoint is absent:

- skip the dependent row;
- add `_bioextract.validation_issue` with
  `issue_code=foreign_key_violation`;
- publish as `passed_with_warnings`.

Partial constructors distinguish an intentionally absent capability from a
broken complete release. Unknown optional entry fields are preserved only when
the schema explicitly owns them; otherwise they are ignored and covered by
parser compatibility tests.

## Implementation Slices

### Slice 1: streaming canonical parser

- add `src/bioextract/kegg/metabolic/`;
- parse list files and the included scalar entity fields;
- parse global relation TSVs;
- parse equations with symbolic coefficients;
- parse module definitions into an AST;
- test multi-record batches, continuation lines, obsolete EC entries, C/G
  participants, and malformed records.

### Slice 2: publication and validation

- add `from_metabolic_release()` and `from_metabolic_files()`;
- publish `kegg-metabolic-v0.1` to one DuckDB;
- create metadata schema v3 and capability metadata;
- verify staging cleanup and old-target protection.

### Slice 3: publication-backed domain access

- add `from_duckdb()` and read-only `connect()`;
- extend `select_ids()`/`select_groups()` for metabolic namespaces;
- implement reaction-centered extraction and unmatched-ID accounting;
- implement exact module block evaluation.

### Slice 4: real-snapshot acceptance

- build:

  ```text
  /cephfs_data/genostack_v3/genostack_php/public_file_data/database/bioinfo/resources/kegg/metabolic/2026-07/tidy/kegg.duckdb
  ```

- verify all source entry and relation counts;
- smoke C number, ChEBI, PubChem, R number, Rhea, EC, KO, module, and pathway
  selections;
- verify grouped isolation and unmatched IDs;
- directly join KEGG ChEBI references to ChEBI and KEGG Rhea references to
  Rhea without string construction or casts;
- benchmark build memory, elapsed time, file size, point selection, and a
  compound-to-pathway traversal.

## Compatibility

The existing BRITE JSON and organism-mapping constructors and Parquet schemas
remain unchanged. Metabolic support is additive under `KEGGDatabase`.

`kegg-metabolic-v0.1` is a separate resource schema. It does not reinterpret
`kegg-brite-tidy-v0.1` or `kegg-mapping-v0.1`, and it does not make their
independent Parquet files obsolete.

KEGG access and redistribution terms remain the caller's responsibility.
`bioextract` consumes caller-supplied local files and performs no network
requests.

## Completion Record

Completed locally on 2026-07-30. The implementation uses exact and
archive-safe release-layout discovery, streaming NDJSON relation spools, the
shared metadata-v3 DuckDB publication lifecycle, capability-to-inventory
validation, validated read-only reopening, reaction-centered PascalCase
extractors, transferred/obsolete EC semantics, and recursive exact module
block evaluation with cycle detection. Focused KEGG tests and `pdm run check`
pass.

The real 2026-07 CephFS layout resolved the expected four lists, four
entry-batch collections, and seven global relations with batch counts
1,962/1,246/835/58. The formal publication was atomically written to
`kegg/metabolic/2026-07/tidy/kegg.duckdb`: 25 canonical tables, 444,929 rows,
17,313,792 bytes, legacy metadata schema v2 (an old artifact, not the current
v3 writer contract), and two persisted
`foreign_key_violation` warnings for absent compound participant `C23109`.
Publication with 4,112 source SHA-256 values took 2:08.57 wall-clock time and
peaked at 1,476,904 KiB RSS.

Acceptance covered all nine identifier namespaces, grouped unmatched-ID
isolation, all-module evaluation, read-only native SQL, and direct shared-ID
joins to the ChEBI and Rhea publications. The reproducible inventory, probes,
and warm-query observations are recorded in the
[2026-07 benchmark](../benchmarks/20260730-kegg-metabolic-2026-07-benchmark.md).
