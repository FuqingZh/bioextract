# bioextract Documentation

Version: v1.0
Date: 2026-08-04
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
| STRINGdb and OmniPath | Root [README](../README.md), source code, and tests | [`STRINGdb integration`](../tests/integration/stringdb/test_database.py), [`OmniPath integration`](../tests/integration/omnipath/test_database.py) | — |

## Plans And History

The completed
[Storage And Publication Convergence Plan](implementation-plan/20260803-v1.0-storage-publication-convergence-implementation-plan.md)
records the accepted direct-access versus materialization boundary, DuckDB-only
bioextract publication contract, metadata v1 reset, former-Parquet resource
migrations, eggNOG SQLite behavior, GO subcell cleanup, and acceptance of the
nine maintained formal publications. All five slices are complete.

The current
[Test Suite Layering Plan](implementation-plan/20260731-v1.0-test-suite-layering.md)
defines the unit, contract, integration, and external-snapshot smoke
boundaries, fixture ownership, migration slices, and unchanged 482-case
baseline.

The completed
[Public API And Grouped Selection Convergence Plan](implementation-plan/20260731-v1.0-public-api-grouped-selection-convergence.md)
records the top-level database import contract, shared unique-ID resolution and
group fan-out model, internal behavior-type boundary, and unchanged
publication schemas.

The current
[Constructor Source Convergence Plan](implementation-plan/20260731-v1.0-constructor-source-convergence.md)
defines the accepted three-phase constructor naming contract and the eggNOG,
STRING, WikiPathways, and Rhea delivery boundaries.

The completed [UniProtKB Domain Access Plan](implementation-plan/20260730-uniprot-knowledgebase-domain-access.md)
records the Swiss-Prot DuckDB boundary and version-contract convergence.

The completed
[KEGG Metabolic Domain Access Plan](implementation-plan/20260730-kegg-metabolic-domain-access.md)
records the compound, reaction, enzyme, module, cross-reference, publication,
and reaction-centered selection boundary.

The completed
[ChEBI Domain Access And Publication Plan](implementation-plan/20260729-chebi-domain-access.md)
records the portable ChEBI publication, exact compound
selection, stable profile extraction, and cycle-safe relation traversal.

The completed
[Domain Contract Convergence Plan](implementation-plan/20260729-v1.0-domain-contract-convergence.md)
records the dependency-ordered public-output, naming, provenance,
limits-removal, and resource-API migration work.

Completed plans with reusable design rationale are retained under
[`archive/`](archive/). Superseded GO term-selection and UniProt subcellular
location plans point to their current architecture authorities. Historical PR
descriptions are not project documentation; Git history remains their archive.
