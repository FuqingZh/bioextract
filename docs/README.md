# bioextract Documentation

Version: v1.0
Date: 2026-07-14
Status: current

This directory contains durable architecture, data contracts, test standards,
and reproducible benchmark baselines. Public usage starts in the repository
[README](../README.md); implementation details in source code remain authoritative
when a domain does not yet have a dedicated architecture document.

## Reading Order

1. Read the shared [Tidy Dataset Contract](architecture/tidy-dataset-contract.md)
   before changing generated Parquet or manifest behavior.
2. Read the resource-specific architecture document before changing a public DB
   handle, schema, selection rule, or write contract.
3. Read the corresponding test standard and, where present, the real-snapshot
   benchmark before changing validation or publication behavior.

## Authority Map

| Area | Architecture or contract | Test standard | Benchmark |
| --- | --- | --- | --- |
| Shared tidy output | [Tidy Dataset Contract](architecture/tidy-dataset-contract.md) | Full pytest suite | — |
| GO ontology and KEGG BRITE | [GO and KEGG Tidy](architecture/go-kegg-tidy.md) | [GO and KEGG Tidy](testing/go-kegg-tidy.md) | — |
| KEGG organism mapping | [KEGG Mapping DB](architecture/kegg-mapping-db.md) | [KEGG Mapping DB](testing/kegg-mapping-db.md) | — |
| Reactome | [ReactomeDb](architecture/reactome-db.md) | [ReactomeDb](testing/reactome-db.md) | — |
| WikiPathways | [WikiPathwaysDb](architecture/wikipathways-db.md) | [WikiPathwaysDb](testing/wikipathways-db.md) | — |
| eggNOG | [EggnogDb](architecture/eggnog-db.md) | [EggnogDb](testing/eggnog-db.md) | [5.0.2 baseline](benchmarks/20260608-v1.0-eggnog-5.0.2-benchmark.md) |
| InterPro and Pfam | [InterProDb](architecture/interpro-db.md) | [InterProDb](testing/interpro-db.md) | [108.0 baseline](benchmarks/20260714-v1.0-interpro-108-benchmark.md) |
| UniProt | [UniprotDb](architecture/uniprot-db.md) | [UniprotDb](testing/uniprot-db.md) | — |
| STRINGdb and OmniPath | Root [README](../README.md), source code, and tests | [`test_stringdb.py`](../tests/test_stringdb.py), [`test_omnipath.py`](../tests/test_omnipath.py) | — |

## Plans And History

There is no active implementation plan in this directory. Create
`docs/implementation-plan/` only when an approved, unfinished plan needs a
durable execution boundary.

Completed plans with reusable design rationale are retained under
[`archive/`](archive/). Superseded GO term-selection and UniProt subcellular
location plans point to their current architecture authorities. Historical PR
descriptions are not project documentation; Git history remains their archive.
