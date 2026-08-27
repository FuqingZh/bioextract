# KEGG Mapping Database Architecture

Version: v2.3
Date: 2026-08-27
Status: current

## Boundary

`KEGGDatabase` reads caller-supplied local KEGG mapping files. It does not
download data, infer releases from paths, read biofetch manifests, or require a
complete KEGG release. A directory containing 200 valid organism directories
is a valid all-available build of those 200 members.

## Construction And Scope

```python
db = KEGGDatabase.from_mapping_directory(
    "/data/kegg/mapping/2026-06/raw",
    release_version="2026-06",
)
hsa = db.with_organisms(["hsa"])

partial = KEGGDatabase.from_mapping_files(
    "/data/kegg/mapping/2026-06/raw/hsa",
    organism_code="hsa",
)
```

`from_mapping_directory(source, *, organism_list, ko_pathway,
release_version)` treats `source` as the direct multi-organism discovery root.
It enumerates only immediate children matching `^[a-z]{3,4}$`; `ko` and
`organism` are reserved global directories. Each discovered organism must
contain the rectangular five-role profile.

`from_mapping_files(source=None, *, organism_code, ...)` treats `source` as one
organism's direct directory. Explicit role files replace conventional
children. Missing roles are unavailable capabilities, and at least one
organism role is required. This factory never climbs to a parent to find global
roles.

`with_organisms()` is the physical-pruning API. A source handle opens only the
selected organism directories. A publication handle may narrow to members in
its `organism` table but cannot expand beyond them.

## Lazy Domain Relations

Every tabular method returns a native replayable `polars.LazyFrame`:

```python
db.organisms()
db.gene_annotations()
db.ko_annotations()
db.gene_pathways()
db.ko_pathways()
db.gene_pathways_via_ko()

selection = db.select_ids(ids, namespace="uniprot")
selection.matches()
selection.gene_annotations()
selection.gene_pathways()
selection.gene_pathways_via_ko()
selection.unmatched_ids()
```

There are no hidden row limits, eager biological caches, `extract_*` aliases,
or `mappings()` compatibility method. Callers use native Polars `collect()`,
`sink_*`, `explode()`, and `unnest()`.

Direct gene-pathway observations and gene-to-KO-to-pathway traversal are
separate relations. The latter retains `ko_id` inside every evidence struct.

## Selection Input Identity

The mapping `uniprot` namespace accepts a trimmed plain accession, its
existing `up:` form, or one complete `sp|accession|entry_name` /
`tr|accession|entry_name` representation. The resulting accession is public
`input_id` and the deduplicated membership key. Any other pipe-bearing
UniProt caller value raises `ValueError`.

`ncbi_gene` retains its existing optional `ncbi-geneid:` rule and otherwise
preserves trimmed text. `kegg_gene` still requires a full organism prefix.
Neither non-UniProt namespace uses the UniProt representation parser, and
official mapping-file rows continue through their resource-owned parsers.

## Aggregate Schemas

Every mapping publication contains exactly:

- `organism(organism_code, genome_id, organism_name, taxonomy_lineage)`;
- `gene_annotation`, one row per `(organism_code, kegg_gene_id)`, with scalar
  gene attributes and nested `uniprot_mappings`, `ncbi_gene_mappings`,
  `ko_mappings`, and `pathway_mappings`; and
- `ko_annotation`, one row per `ko_id` with nested pathway mappings.

Nested relations are `List[Struct]`. `null` means the required source role was
unavailable, `[]` means it was available but no observation existed, and a
non-empty list contains distinct locally sorted observations. Top-level row
order is not a caller contract.

The gene universe is the union of gene IDs in every available organism role.
The KO universe is the union of organism gene-to-KO observations and the
optional global KO-pathway relation.

## Parsing And Integrity

Files are UTF-8, headerless TSV. Blank lines are ignored; non-blank rows require
the exact role-specific column count. Fatal errors identify the logical role,
path when available, and a cause category; line details are diagnostic rather
than a scheduling or compatibility promise. Gene IDs must match the current
organism prefix, KO and pathway IDs use their fixed KEGG grammars, and NCBI
Gene values remain strings. Exact duplicate relationship rows are idempotent
and removed silently; conflicting metadata remains fatal.

`gene_list.tsv` has four fields: gene ID, gene type, opaque genomic position,
and display text. Display text is split once at the first semicolon; comma-
separated names become one symbol plus sorted aliases.

## Publication

`write_duckdb()` uses one DuckDB-native SQL coordinator. It enumerates the build
scope once, records a declarative input inventory, scans global roles once, and
processes organism roles in private hash windows. Each window bulk-scans an
explicit path list, validates and pre-aggregates its many-side relations, then
releases the window state before the next one. The writer uses an adjacent
staging database and per-publication temporary spill directory, with internal
thread, memory, and temporary-size bounds; these are execution safeguards, not
row limits or public tuning parameters. It embeds metadata, validates read-only,
drops all private relations, and atomically replaces the destination. It does
not run the three public lazy relations independently or create per-organism
Parquet publication files.

Identity is fixed:

```text
bioextract.metadata_schema_version = "2"
bioextract.resource_schema_version = "kegg-mapping-v1.0"
bioextract.source_schema_profile = "kegg-organism-mapping-files-v2"
bioextract.scope = "mapping"
```

Seven exact boolean capability keys describe the five organism roles and two
global roles. `_bioextract.source_file` is the sole declarative source
inventory. `bytes` and `sha256` are nullable; the writer does not add a
provenance-only stat or content pass. No `manifest.lock` or sidecar is read or
written.

Old wide-table mapping publications and metadata v1 are intentionally rejected
without aliases or migration views.

A reopened mapping handle pins the filesystem identity that remained stable
across profile and publication validation. Every native `connect()`, selection
match, and replayable publication relation checks that identity before opening
DuckDB and after a successful read. `with_organisms()` retains the same pinned
identity. If the path is removed or atomically replaced, the old handle fails
with `IntegrityError` and the caller must explicitly reopen it with
`from_duckdb()`.
