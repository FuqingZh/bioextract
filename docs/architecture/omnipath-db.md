# OmniPath Database Architecture

`OmniPathDatabase` is a direct, local-snapshot adapter. It accepts one or both
official relation files and does not download data or create a separate
publication. `available_resources` reports the roles present in the supplied
snapshot; `scan_enzsub()` and `scan_interactions()` expose the official columns
and order lazily.

The current source profile represents protein identifiers, so selections use
`namespace="protein"`. The explicit namespace parameter leaves room for a
future source profile without guessing from filenames or identifier values.

Enzyme-substrate selections preserve the official `enzyme`, `substrate`,
`residue_type`, `residue_offset`, and `modification` columns. `target_site` is
a derived convenience field and never replaces its source fields. Interaction
selections preserve the official endpoint and nullable evidence columns. Group
selection resolves globally unique identifiers once, then fans matched rows
back through membership; records differing in evidence are not collapsed.

```python
from bioextract import OmniPathDatabase

database = OmniPathDatabase.from_files(
    enzsub="enzsub.tsv.gz",
    interactions="interactions.tsv.gz",
)
selection = database.select_groups(
    {"case": ["P31749"], "control": ["MAPK1"]},
    namespace="protein",
)
lf_enzsub = selection.enzsub()
lf_interactions = selection.with_interactions().interactions()

# Native Polars execution remains the caller's choice.
df_enzsub = lf_enzsub.collect()
```

Convenience existence probes are intentionally not part of the public handle;
callers can inspect the lazy source scans or use a selection. This keeps source
inspection and domain extraction on one composable Polars boundary.
