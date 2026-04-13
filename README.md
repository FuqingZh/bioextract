# bioextract

Library-first extraction helpers for bioinformatics resource snapshots.

## Install

- `pip install bioextract`

## STRINGdb

```python
from bioextract.stringdb import StringDb, StringResourceLimits

selection = (
    StringDb.from_files(
        file_aliases="9606.protein.aliases.v12.0.txt.gz",
        file_links="9606.protein.links.v12.0.txt.gz",
        limits=StringResourceLimits(num_input_ids_max=50_000),
    )
    .select_ids(["P04637", "EGFR", "CDK2"])
    .with_score_min(400)
)

df_mapping = selection.extract_string_mapping()
df_unmapped = selection.extract_unmapped_input_ids()
df_edges = selection.extract_edges()

print(df_mapping)
print(df_unmapped)
print(df_edges)
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
