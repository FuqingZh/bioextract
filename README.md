# bioextract

Library-first extraction helpers for bioinformatics resource snapshots.

## Install

- `pip install bioextract`

## STRINGdb

```python
from bioextract.stringdb import StringDb, StringResourceLimits

result = (
    StringDb.from_files(
        file_aliases="9606.protein.aliases.v12.0.txt.gz",
        file_links="9606.protein.links.v12.0.txt.gz",
        limits=StringResourceLimits(num_input_ids_max=50_000),
    )
    .select_ids(["P04637", "EGFR", "CDK2"])
    .with_score_min(400)
    .extract_result()
)

print(result.df_edges)
print(result.metrics.num_edges)
```

```python
from bioextract.stringdb import StringDb

df_group_edges = (
    StringDb.from_files(
        file_aliases="9606.protein.aliases.v12.0.txt.gz",
        file_links="9606.protein.links.v12.0.txt.gz",
    )
    .select_groups(
        {
            "TumorA": ["TP53", "EGFR"],
            "TumorB": ["CDK2", "TP53"],
        }
    )
    .with_score_min(400)
    .extract_edges()
)
```

## Development

- `PYTHONPATH=src pytest`
- `PYTHONPATH=src python scripts/benchmark_stringdb.py`
