# UniProt idmapping scan benchmark

Version: v1.0
Date: 2026-08-12
Status: bounded fixture decision record

## Scope

This is a bounded parser/execution check for the existing
`UniProtDatabase.scan_mapping(taxon_ids=...)` path. It is not a production
snapshot benchmark: no selected-ID query workload has been supplied with a
publication identity and caller SLA, so the result is not used to invent a
second selection API.

## Reproduction

The check creates a temporary gzip idmapping file with 100,000 rows and runs:

```python
database.scan_mapping(taxon_ids=None).collect()
database.scan_mapping(taxon_ids=["9606"]).collect()
database.scan_mapping(taxon_ids=["10090"]).collect()
```

The command was run with `pdm run python` from the repository working tree on
2026-08-12. The fixture has one taxon (`9606`) and the standard 22 selected
columns. RSS is sampled from `/proc/self/status` immediately before and after
each collection; it is an observation, not a memory limit.

## Observed result

| query | rows | elapsed | RSS delta |
| --- | ---: | ---: | ---: |
| no taxon filter | 100,000 | 0.158 s | 143.5 MiB |
| `taxon_ids=["9606"]` | 100,000 | 0.051 s | 43.3 MiB |
| `taxon_ids=["10090"]` | 0 | 0.034 s | 39.2 MiB |

The taxon predicate is effective on the bounded fixture, and every execution
returns a replayable native `pl.LazyFrame` with fresh source ownership.

## Decision

Keep `scan_mapping(taxon_ids=...)` as the physical idmapping API. Do not add a
publication-backed selected-ID API in this release. Revisit only when a real
caller workload records snapshot identity, selected-ID cardinality, cold/warm
timings, peak RSS, checksum-equivalent row counts, and an explicit operational
budget showing that the existing scan misses it.

The bounded fixture does not prove production-scale behavior; that residual
risk is recorded rather than hidden behind a new API or an eager cache.
