# bioextract Documentation

Version: v1.0
Date: 2026-08-05
Status: current

This directory contains durable architecture, data contracts, test standards,
and reproducible benchmark baselines. Public usage starts in the repository
[README](../README.md); implementation details in source code remain authoritative
when a domain does not yet have a dedicated architecture document.

## Reading Order

1. Read the repository-wide
   [Domain Access Architecture](architecture/20260729-v1.0-domain-access-architecture.md)
   before adding a resource, public query method, or storage strategy.
2. Read the shared [Materialized Dataset Contract](architecture/tidy-dataset-contract.md)
   before changing DuckDB publication or embedded provenance behavior.
3. Read the resource-specific architecture document before changing a public DB
   handle, schema, selection rule, or write contract.
4. Read the repository-wide [Test Standard](testing/README.md), then the
   corresponding resource test standard and, where present, the real-snapshot
   benchmark before changing validation or publication behavior.
5. Read the repository [AO and delivery norms](runbooks/ao-delivery.md)
   before starting PR-bound work or changing repository gates.

## Directory Map

- `architecture/` — current system boundaries and materialized-data contracts.
- `implementation-plans/` — accepted execution plans and their historical
  closeout records.
- `testing/` — reusable validation standards and resource-specific checks.
- `benchmarks/` — reproducible measurement baselines.
- `runbooks/` — controlled delivery and operational procedures.
- `archive/` — superseded plans retained for traceability.

## Authority Map

| Area | Architecture or contract | Test standard | Benchmark |
| --- | --- | --- | --- |
| Repository purpose and boundaries | [Domain Access Architecture](architecture/20260729-v1.0-domain-access-architecture.md) | [Repository test standard](testing/README.md) | — |
| Shared publication output | [Materialized Dataset Contract](architecture/tidy-dataset-contract.md) | [Repository test standard](testing/README.md) | — |
| ChEBI and ChemOnt | [ChEBIDatabase](architecture/chebi-db.md) | [ChEBIDatabase](testing/chebi-db.md) | — |
| GO ontology and KEGG BRITE | [GO and KEGG Tidy](architecture/go-kegg-tidy.md) | [GO and KEGG Tidy](testing/go-kegg-tidy.md) | — |
| KEGG organism mapping | [KEGG Mapping DB](architecture/kegg-mapping-db.md) | [KEGG Mapping DB](testing/kegg-mapping-db.md) | — |
| KEGG metabolic | [KEGG Metabolic Database](architecture/kegg-metabolic-db.md) | [KEGG Metabolic Database](testing/kegg-metabolic-db.md) | [2026-07 baseline](benchmarks/20260730-kegg-metabolic-2026-07-benchmark.md) |
| Reactome | [ReactomeDatabase](architecture/reactome-db.md) | [ReactomeDatabase](testing/reactome-db.md) | — |
| WikiPathways | [WikiPathwaysDatabase](architecture/wikipathways-db.md) | [WikiPathwaysDatabase](testing/wikipathways-db.md) | — |
| eggNOG | [EggNOGDatabase](architecture/eggnog-db.md) | [EggNOGDatabase](testing/eggnog-db.md) | [5.0.2 baseline](benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md) |
| InterPro and Pfam | [InterProDatabase](architecture/interpro-db.md) | [InterProDatabase](testing/interpro-db.md) | [108.0 baseline](benchmarks/20260714-v1.0-interpro-108-benchmark.md) |
| UniProt | [UniProtDatabase](architecture/uniprot-db.md) | [UniProtDatabase](testing/uniprot-db.md) | [2026_01 baseline](benchmarks/20260731-uniprot-kb-2026_01-benchmark.md) |
| Rhea | [RheaDatabase](architecture/rhea-db.md) | [RheaDatabase](testing/rhea-db.md) | — |
| STRINGdb | Root [README](../README.md), source code, and tests | [`STRINGdb integration`](../tests/integration/stringdb/test_database.py) | — |
| OmniPath | [OmniPath architecture](architecture/omnipath-db.md) | [OmniPath test standard](testing/omnipath-db.md) | — |

## Plans And History

The active
[Runtime And Constructor Convergence Plan](implementation-plans/20260805-v1.0-runtime-and-constructor-convergence.md)
defines the repository validation thread-budget boundary, the separate
host/AO resource-control ownership, the common `source` overlay contract, and
the dependency-ordered Rhea, KEGG metabolic, and ChEBI constructor migrations.

The active
[Column Lineage And Public Schema Convergence Plan](implementation-plans/20260804-v1.0-column-lineage-public-schema-convergence.md)
defines the official-header versus derived-column ownership rule, shared
`snake_case` selection protocol, `from_*` lifecycle naming, resource-specific
migration slices, and the publication rebuild and verification boundary.

The completed
[Publication Inspection Plan](implementation-plans/20260804-v1.0-publication-inspection-implementation-plan.md)
records the narrow, read-only API for inspecting exactly one caller-supplied
local bioextract DuckDB publication, its metadata-v1 validation and immutable
result contract, and the boundary that leaves optional discovery and catalogs
to external callers.

The completed
[Storage And Publication Convergence Plan](implementation-plans/20260803-v1.0-storage-publication-convergence-implementation-plan.md)
records the accepted direct-access versus materialization boundary, DuckDB-only
bioextract publication contract, metadata v1 reset, former-Parquet resource
migrations, eggNOG SQLite behavior, GO subcell cleanup, and acceptance of the
nine maintained formal publications. All five slices are complete.

The completed
[Test Suite Layering Plan](implementation-plans/20260731-v1.0-test-suite-layering.md)
records the unit, contract, integration, and external-snapshot smoke
boundaries, fixture ownership, migration slices, and the growth from its
482-case pre-migration baseline to the current 598-case hermetic gate.

The completed
[Public API And Grouped Selection Convergence Plan](implementation-plans/20260731-v1.0-public-api-grouped-selection-convergence.md)
records the top-level database import contract, shared unique-ID resolution and
group fan-out model, internal behavior-type boundary, and unchanged
publication schemas.

The superseded, completed
[Constructor Source Convergence Plan](archive/20260731-v1.0-constructor-source-convergence.md)
records the implemented eggNOG, STRING, WikiPathways, and first Rhea
constructor delivery boundaries. The active runtime and constructor plan
replaces its frozen Rhea, KEGG metabolic, and ChEBI follow-on decisions.

The completed [UniProtKB Domain Access Plan](implementation-plans/20260730-uniprot-knowledgebase-domain-access.md)
records the Swiss-Prot DuckDB boundary and version-contract convergence.

The completed
[KEGG Metabolic Domain Access Plan](implementation-plans/20260730-kegg-metabolic-domain-access.md)
records the compound, reaction, enzyme, module, cross-reference, publication,
and reaction-centered selection boundary.

The completed
[ChEBI Domain Access And Publication Plan](implementation-plans/20260729-chebi-domain-access.md)
records the portable ChEBI publication, exact compound
selection, stable profile extraction, and cycle-safe relation traversal.

The completed
[Domain Contract Convergence Plan](implementation-plans/20260729-v1.0-domain-contract-convergence.md)
records the dependency-ordered public-output, naming, provenance,
limits-removal, and resource-API migration work.

Superseded plans with reusable design rationale are retained under
[`archive/`](archive/). The archived constructor, GO term-selection, and
UniProt subcellular-location plans point to their current or follow-up
authorities. Historical PR descriptions are not project documentation; Git
history remains their archive.
