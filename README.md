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

## OmniPath

```python
from bioextract.omnipath import OmniPathDb

selection = (
    OmniPathDb.from_files(
        file_enzsub="enzsub.tsv.gz",
        file_interactions="interactions.tsv.gz",
    )
    .select_ids(["P31749", "AKT1", "BAD"])
    .with_enzsub()
)

df_enzsub = selection.extract_enzsub()
df_unmapped = selection.extract_unmapped_input_ids()

print(df_enzsub)
print(df_unmapped)
```

```python
from bioextract.omnipath import OmniPathDb

df_group_interactions = (
    OmniPathDb.from_files(file_interactions="interactions.tsv.gz")
    .select_groups(
        {
            "TumorA": ["AKT1", "MTOR"],
            "TumorB": ["EGFR", "ERBB2"],
        }
    )
    .with_interactions()
    .extract_interactions()
)
```

## Development

- `PYTHONPATH=src pytest`
- `PYTHONPATH=src python scripts/benchmark_stringdb.py`

## Release

- GitHub Actions now provides:
  - `.github/workflows/py-ci.yml` for test-and-build checks on push and pull request
  - `.github/workflows/publish.yml` for tag-triggered PyPI publishing
- Release tags must be canonical PEP 440 versions such as `0.1.1`
- The publish workflow expects PyPI trusted publishing to be configured for the `pypi` environment
