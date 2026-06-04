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

## GO

```python
from bioextract.go import GoDb

tidy = GoDb.from_obo("go-basic.obo").build_tidy()

df_term = tidy.frames["term"]
df_edge = tidy.frames["edge"]
df_ancestor = tidy.frames["ancestor_all"]
df_subcell = GoDb.from_obo("go-basic.obo").extract_subcell()

report = tidy.write("out/go-basic")
```

`GoDb.from_obo(...).write_tidy("out/go-basic")` is also available as a
convenience wrapper when only persisted parquet outputs are needed.
Pass `should_write_manifest=True` to also write `manifest.json`.
`GoDb.from_obo(...).write_subcell("out/subcell.parquet")` writes non-obsolete
cellular component terms as a subcellular-location table.

## KEGG

```python
from bioextract.kegg import KeggDb

tidy = KeggDb.from_brite_json("br08901.json").build_tidy()

df_pathway = tidy.frames["pathway"]

report = tidy.write("out/br08901")
```

The GO and KEGG tidy writers emit flat parquet files by default. See
`docs/architecture/go-kegg-tidy.md`.

## Reactome

```python
from bioextract.reactome import ReactomeDb

db = ReactomeDb.from_files(
    file_uniprot2reactome="UniProt2Reactome.txt",
    file_pathways="ReactomePathways.txt",
    file_relations="ReactomePathwaysRelation.txt",
)

selection = db.with_species("Homo sapiens").select_ids(["P04637", "Q9Y243"])

df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()
df_term2gene = db.with_species("Homo sapiens").extract_term2gene()
df_term2name = db.with_species("Homo sapiens").extract_term2name()
```

`ReactomeDb` reads local Reactome mapping files and emits annotation tables plus
standard enrichment inputs. The three raw files are composable: mapping-only
snapshots can still emit `mapping` and `term2gene`, pathways-only snapshots can
emit `pathway` and `term2name`, and relation extraction uses the relation file.
It does not call Reactome web services or calculate enrichment p-values.

## WikiPathways

```python
from bioextract.wikipathways import WikiPathwaysDb

db = WikiPathwaysDb.from_gmt(
    "wikipathways-20260510-gmt-Homo_sapiens.gmt",
    species="Homo sapiens",
)

df_pathway = db.extract_pathway()
df_term2gene = db.extract_term2gene()
df_term2name = db.extract_term2name()

selection = db.select_ids(["2687", "435", "MISSING"])
df_mapping = selection.extract_mapping()
df_unmapped = selection.extract_unmapped_input_ids()

report = db.write_tidy("out/wikipathways-hsa", should_write_manifest=True)
```

`WikiPathwaysDb` reads local WikiPathways GMT files. GMT gene content is treated
as NCBI Entrez Gene IDs; the library does not perform identifier conversion or
calculate enrichment p-values.

## UniProt

```python
from bioextract.uniprot import UniprotDb

db = UniprotDb.from_files(
    file_idmapping_selected="idmapping_selected.tab.gz",
)

df_hsa = db.with_taxids("9606").extract_mapping()

report = db.with_taxids("9606", "10090").write_tidy(
    "out/uniprot-idmapping",
    should_write_manifest=True,
)
```

`UniprotDb` reads raw UniProt `idmapping_selected.tab(.gz)`, single parquet
files, or hive parquet dataset directories. Tidy writing emits a canonical
`mapping.parquet`; all-taxid export requires `should_allow_all=True`.
Use `policy_existing="overwrite"` or `policy_existing="skip"` when the output
directory already exists.

## Development

- `PYTHONPATH=src pytest`
- `PYTHONPATH=src python scripts/benchmark_stringdb.py`

## Release

- GitHub Actions now provides:
  - `.github/workflows/py-ci.yml` for test-and-build checks on push and pull request
  - `.github/workflows/publish.yml` for tag-triggered PyPI publishing
- Release tags must be canonical PEP 440 versions such as `0.1.1`
- The publish workflow expects PyPI trusted publishing to be configured for the `pypi` environment
