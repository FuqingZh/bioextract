# KEGG Metabolic Database Test Standard

Version: v0.1
Date: 2026-08-07
Status: current

The focused suites are `tests/unit/kegg/test_metabolic_parsers.py`,
`tests/contract/resources/kegg/test_metabolic_publication_contract.py`, and
`tests/integration/kegg/test_metabolic.py`. They use temporary local fixtures
only and verify:

- streaming multi-record compound, reaction, enzyme, and module parsing;
- partial explicit profiles and `source` discovery from a release root, `raw`,
  nested `raw`, zip, and tar;
- zero/one/multiple layout handling, whole-role overlays, required-role
  completeness, empty-fileset rejection, deterministic entry-batch order, and
  duplicate physical-file rejection;
- optional-list discovery and conditional list/entry parity;
- final directory provenance and archive-plus-overlay provenance;
- symbolic/numeric coefficients and C/G equation participants;
- polymer/side suffix qualifiers without corrupting stoichiometric coefficients;
- active and real `ENTRY EC ... Obsolete Enzyme` records, recursive replacement
  chains, deleted entries, invalid targets, and exact obsolete selection;
- module `--` placeholders, M-number references, reference cycles, and blank
  reaction separators;
- normalized ChEBI and Rhea cross-references;
- all canonical global relation roles and metadata schema v1;
- validated read-only reopening and rejection of wrong resource identity;
- reaction-centered single and grouped selections, extraction terminals, and
  unmatched-ID accounting;
- partial-publication capability failures and namespace availability based on
  actual cross-reference values;
- exact module block evaluation;
- preservation of existing KEGG BRITE and organism-mapping behavior.

Run the focused standard with:

```console
pdm run pytest tests/unit/kegg tests/contract/resources/kegg \
  tests/integration/kegg
```

Before handoff, run the repository-wide non-mutating gate:

```console
pdm run check
```

The accepted real-snapshot inventory, cross-resource joins, domain probes, and
performance observations are recorded in the
[KEGG metabolic 2026-07 baseline](../benchmarks/20260730-kegg-metabolic-2026-07-benchmark.md).
