# bioextract Documentation

Version: v1.8
Date: 2026-08-26
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
5. Read the repository [delivery boundaries](runbooks/repository-delivery.md)
   before starting PR-bound work or changing repository gates.

The [2026-08-26 tree review](architecture/20260826-v1.0-bioextract-tree-review.md)
is observational. It does not replace the domain architecture or tidy-dataset
contract. Remediation order is the
[tree-review remediation plan](implementation-plans/20260826-v1.0-tree-review-remediation-implementation-plan.md).

## Recorded Tree Review

The tree review of commit `64a8802` records that the independence boundary and
metadata-v2 publication path hold in executable code, while the root README
still teaches metadata v1 and deleted `extract_*` terminals, and package
identity is split across `pyproject.toml` `0.1.0`, `dist/` `0.5.0`, and later
plan versions. It does not authorize CephFS publication or a new resource.

## Directory Map

- `architecture/` — current system boundaries and materialized-data contracts.
  Recorded tree reviews in this directory are observational and are not current
  contracts.
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
| Cross-resource lazy relations and execution performance | [Lazy Domain Query Convergence Plan](implementation-plans/20260812-v1.0-lazy-domain-query-convergence-implementation-plan.md); [Execution and Publication Performance Convergence Plan](implementation-plans/20260814-v1.0-cross-resource-execution-performance-convergence-implementation-plan.md) | [Repository test standard](testing/README.md) | [joint_omics 23362 baseline](benchmarks/20260813-v1.0-joint-omics-23362-lazy-relation-baseline.md); [P2](benchmarks/20260814-v1.0-cross-resource-p2-candidates.md), [P3](benchmarks/20260814-v1.0-cross-resource-p3-selected-lookup.md), [P4](benchmarks/20260814-v1.0-cross-resource-p4-publication-loader.md), [v2-compatible joint_omics experiment](benchmarks/20260814-v1.0-cross-resource-p5-v2-compatible-joint-omics.md) |
| Shared publication output | [Materialized Dataset Contract](architecture/tidy-dataset-contract.md); [metadata v2 declarative-provenance amendment](implementation-plans/20260813-v1.0-publication-metadata-v2-migration-plan.md) | [Repository test standard](testing/README.md) | — |
| Formal metadata-v2 publications and release catalog | [20260817-r1 Formal Publication Release Plan](implementation-plans/20260817-v1.0-formal-publication-release-implementation-plan.md) | [Repository test standard](testing/README.md); resource-specific standards | [bioextract-owned reference workload](benchmarks/20260817-v1.0-formal-release-reference-workload.toml); [candidate evidence](benchmarks/20260817-v1.0-formal-release-candidate-evidence.toml); [A-only noise](benchmarks/20260817-v1.0-formal-release-reference-a-noise.json); [interleaved A/B](benchmarks/20260817-v1.0-formal-release-reference-ab.json) |
| ChEBI and ChemOnt | [ChEBIDatabase](architecture/chebi-db.md) | [ChEBIDatabase](testing/chebi-db.md) | — |
| GO ontology and KEGG BRITE | [GO and KEGG Tidy](architecture/go-kegg-tidy.md) | [GO and KEGG Tidy](testing/go-kegg-tidy.md) | — |
| KEGG organism mapping | [KEGG Mapping DB](architecture/kegg-mapping-db.md); P3 DuckDB-native publisher accepted on the fixed 1,000-organism scope, with [release held](implementation-plans/20260813-v1.0-kegg-mapping-aggregate-publication-implementation-plan.md) pending its separate release gate | [KEGG Mapping DB](testing/kegg-mapping-db.md) | [2026-06 candidate benchmark](benchmarks/20260813-v1.0-kegg-mapping-2026-06-candidate.md); [full-build evidence](benchmarks/20260817-v1.0-formal-kegg-mapping-build.json); [rno scoped readback](benchmarks/20260817-v1.0-formal-kegg-mapping-rno-validation.json) |
| KEGG metabolic | [KEGG Metabolic Database](architecture/kegg-metabolic-db.md) | [KEGG Metabolic Database](testing/kegg-metabolic-db.md) | [2026-07 baseline](benchmarks/20260730-kegg-metabolic-2026-07-benchmark.md) |
| Reactome | [ReactomeDatabase](architecture/reactome-db.md); [program index](implementation-plans/20260817-v1.0-reactome-mapping-capability-expansion-implementation-plan.md); [P2 mapping matrix](implementation-plans/20260818-v1.0-reactome-mapping-matrix-expansion-implementation-plan.md); [P3 human entities](implementation-plans/20260818-v1.0-reactome-human-entity-pathway-implementation-plan.md); [P4 human GMT](implementation-plans/20260818-v1.0-reactome-human-gene-set-implementation-plan.md) | [ReactomeDatabase](testing/reactome-db.md) | [bounded v96 P1 smoke](benchmarks/20260817-v1.0-reactome-p1-v96-smoke.json); [complete v0.5 smoke](benchmarks/20260818-v1.0-reactome-v05-v96-smoke.json) |
| WikiPathways | [WikiPathwaysDatabase](architecture/wikipathways-db.md) | [WikiPathwaysDatabase](testing/wikipathways-db.md) | — |
| eggNOG | [EggNOGDatabase](architecture/eggnog-db.md) | [EggNOGDatabase](testing/eggnog-db.md) | [5.0.2 baseline](benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md) |
| InterPro and Pfam | [InterProDatabase](architecture/interpro-db.md) | [InterProDatabase](testing/interpro-db.md) | [108.0 baseline](benchmarks/20260714-v1.0-interpro-108-benchmark.md) |
| UniProt | [UniProtDatabase](architecture/uniprot-db.md) | [UniProtDatabase](testing/uniprot-db.md) | [KB 2026_01 baseline](benchmarks/20260731-uniprot-kb-2026_01-benchmark.md); [idmapping scan decision](benchmarks/20260812-uniprot-idmapping-scan-benchmark.md) |
| Rhea | [RheaDatabase](architecture/rhea-db.md) | [RheaDatabase](testing/rhea-db.md) | — |
| STRINGdb | Root [README](../README.md), source code, and tests | [`STRINGdb integration`](../tests/integration/stringdb/test_database.py) | — |
| OmniPath | [OmniPath architecture](architecture/omnipath-db.md) | [OmniPath test standard](testing/omnipath-db.md) | — |

## Plans And History

The active
[tree-review remediation plan](implementation-plans/20260826-v1.0-tree-review-remediation-implementation-plan.md)
orders README/metadata-v2 honesty, package-identity alignment, UniProt ID
normalization, bounded directory discovery, and related integrity fixes. It is
not runtime authority and does not publish formal artifacts.

The accepted
[Reactome Mapping Capability Expansion Plan](implementation-plans/20260817-v1.0-reactome-mapping-capability-expansion-implementation-plan.md)
is now the Reactome program index. P1 delivered explicit
`UniProt2Reactome_All_Levels.txt` source, query, selection, v0.2 publication,
reopen, package 0.6.0, and authorized formal-v96/catalog support. The remaining
implementation authority was ordered as [P2: the ten-file mapping-matrix
expansion](implementation-plans/20260818-v1.0-reactome-mapping-matrix-expansion-implementation-plan.md),
[P3: human Complex/EWAS relations](implementation-plans/20260818-v1.0-reactome-human-entity-pathway-implementation-plan.md),
and [P4: human GMT gene sets](implementation-plans/20260818-v1.0-reactome-human-gene-set-implementation-plan.md),
targeting v0.3, v0.4, and v0.5 respectively. The implementation and temporary
v96 publication smoke now cover the complete v0.5 API; existing UniProt
pathway defaults remain unchanged. The coordinated package, formal Reactome
v96/KEGG publication, and catalog closeout is the remaining release boundary;
downstream activation remains separate.

The accepted
[Formal Metadata v2 And KEGG Mapping Release Plan](implementation-plans/20260817-v1.0-formal-publication-release-implementation-plan.md)
owns the original `20260817-r1` formal resource-delivery baseline: five frozen
bioextract core publications, the then-current Reactome profile under
bioextract 0.5.0, one all-available-organism KEGG mapping DuckDB, and the
immutable release catalog. The earlier Reactome v0.2 and bioextract 0.6.0
formal closeout is recorded in the Reactome program index; the current
coordinated task extends the same release boundary to v0.5. `rno` remains a
reopened-publication selection gate rather than a separate formal artifact,
downstream activation is outside the plan, and the deprecated global
`meta/catalog.toml` is not updated.

The completed
[Publication Metadata v2 Source Inventory Migration Plan](implementation-plans/20260813-v1.0-publication-metadata-v2-migration-plan.md)
removes the duplicated `bioextract.sources` JSON inventory, makes
`_bioextract.source_file` authoritative, treats it as a declarative resolved
input inventory, makes bytes and pre-existing hashes optional without an extra
path/content pass, and requires a coordinated rebuild rather than dual-reading
metadata v1. The nullable-byte and no-hash-option correction is included in the
same pre-release contract slice; no compatibility reader or alias is retained.

The accepted
[KEGG Mapping Aggregate Publication Plan](implementation-plans/20260813-v1.0-kegg-mapping-aggregate-publication-implementation-plan.md)
depends on metadata v2 and replaces the current one-organism wide `mapping`
shape with bounded multi-organism discovery, explicit organism scoping,
aggregate `List(Struct)` gene/KO relations, source/publication lazy parity, and
an evidence gate before any reverse-lookup acceleration table. Its P3
DuckDB-native publisher is accepted on the fixed 100/1,000-organism gates;
package release, a later full build, and formal publication replacement remain
separately authorized.

The completed
[Lazy Domain Query Convergence Plan](implementation-plans/20260812-v1.0-lazy-domain-query-convergence-implementation-plan.md)
defines the repository-wide native `LazyFrame` relation contract, replayable
DuckDB execution ownership, non-Cartesian Rhea `List(Struct)` neighborhoods,
species-safe WikiPathways access, public InterPro Pfam annotations, and the
benchmark gate that precedes any UniProt idmapping selected-ID API. Its final
repository gate passed with 781 tests. The production-scale
[joint_omics 23362 baseline](benchmarks/20260813-v1.0-joint-omics-23362-lazy-relation-baseline.md)
now supplies that workload evidence, confirms shared-execution and physical
pushdown needs, and closes the predecessor's evidence gap without rewriting its
historical decision.

The accepted
[Execution and Publication Performance Convergence Plan](implementation-plans/20260814-v1.0-cross-resource-execution-performance-convergence-implementation-plan.md)
is the follow-up implementation authority. It defines one private
request-aware scan contract with native Polars, DuckDB-to-Arrow,
parser-to-Arrow, and explicit-materialization strategies; compact shared
selection anchors; UniProt idmapping selected access; publication loader
convergence; and executable diagnostic, behavior, benchmark, and gate
self-protection checks. Its P0-P4 implementation slices and bounded P5
repository gate are complete. The P5a Rhea-anchor follow-up on a temporary
v2-compatible joint_omics tree records lower wall time, CPU, and RSS with
stable biological outputs; formal publication and downstream contract work
remain separately authorized. No package or formal publication release is
implied.

The implemented
[GO Ancestor Projection And Publication Query Pushdown Plan](implementation-plans/20260810-v1.0-go-ancestor-query-implementation-plan.md)
records the additive GO term-to-ancestor and GO subset projection API, the
parameterized DuckDB query-pushdown path for publication-backed GO domain
reads, source/publication parity, unchanged nine-table publication contract,
and the boundary that leaves protein membership and enrichment analysis to
downstream applications.

The completed
[Runtime And Constructor Convergence Plan](implementation-plans/20260805-v1.0-runtime-and-constructor-convergence.md)
records the repository validation thread-budget boundary, the separate host/AO
resource-control ownership, gzip FASTA error normalization, the common
`source` overlay contract, and the completed Rhea, KEGG metabolic, and ChEBI
constructor migrations. The implementation slices, independent harness
cleanup, and documentation closeout are merged; the separately authorized
`0.4.0` package-release handoff remains outside the plan.

The completed
[Column Lineage And Public Schema Convergence Plan](implementation-plans/20260804-v1.0-column-lineage-public-schema-convergence.md)
records the official-header versus derived-column ownership rule, shared
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
482-case pre-migration baseline to its 598-case post-migration baseline.

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
