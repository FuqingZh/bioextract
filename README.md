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

go = GoDb.from_obo("go-basic.obo")
tidy = go.build_tidy()

df_term = tidy.frames["term"]
df_edge = tidy.frames["edge"]
df_ancestor = tidy.frames["ancestor_all"]
df_subsets = go.list_subsets()
df_goslim_generic = go.select_terms(subset_id="goslim_generic")
df_subcell = go.extract_subcell()

report = tidy.write("out/go-basic")
```

`GoDb.from_obo(...).write_tidy("out/go-basic")` is also available as a
convenience wrapper when only persisted parquet outputs are needed.
Pass `should_write_manifest=True` to also write `manifest.json`.
`GoDb.from_obo(...).write_subcell("out/subcell.parquet")` writes non-obsolete
cellular component terms as a subcellular-location table.
GO OBO subset memberships are available through `select_terms(subset_id=...)`
and the `subset_membership` tidy frame.

## KEGG

```python
from bioextract.kegg import KeggDb

tidy = KeggDb.from_brite_json("br08901.json").build_tidy()

df_pathway = tidy.frames["pathway"]

report = tidy.write("out/br08901")
```

The GO and KEGG tidy writers emit flat parquet files by default. See the
[GO and KEGG Tidy Architecture](docs/architecture/go-kegg-tidy.md).

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

## eggNOG

```python
from bioextract.eggnog import EggnogDb

db = EggnogDb.from_files(
    file_eggnog_db="eggnog.db.gz",
    file_cog_fun="cog-24.fun.tab",
)

df_mapping = db.extract_mapping()
report = db.write_tidy("out/eggnog", should_write_manifest=True)
```

`EggnogDb` reads local eggNOG mapper SQLite snapshots and optional COG
function lookup tables. It emits a wide protein-to-COG mapping table for
annotation and enrichment-input preparation.

## InterPro

```python
from bioextract.interpro import InterProDb

db = InterProDb.from_mapping_files(
    file_protein2ipr="108.0/raw/protein2ipr.dat.gz",
    file_interpro_xml="108.0/raw/interpro.xml.gz",
)

df_mapping = db.extract_mapping()
report = db.write_tidy("out/interpro", should_write_manifest=True)

pfam_report = db.write_tidy(
    "out/interpro-pfam",
    config="pfam",
    should_write_manifest=True,
)
```

`InterProDb` reads local `protein2ipr` mapping files and optional InterPro XML
metadata, then emits one canonical UniProt-to-InterPro mapping parquet.
`config="pfam"` derives compact UniProt-to-Pfam assets directly from the same
raw snapshot; it does not require the full `mapping.parquet` to be generated
first.

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

For UniProt knowledge-base flat files, `write_eggnog_xref_tidy()` can emit a
canonical UniProt-to-eggNOG xref parquet from `uniprot_sprot.dat(.gz)`.
The same flat-file handle can also extract curated Swiss-Prot subcellular
location comments:

```python
from bioextract.uniprot import UniprotDb

db = UniprotDb.from_dat(
    file_dat="uniprot_sprot.dat.gz",
    source_db="sprot",
)

df_subcell = db.extract_subcellular_location()

report = db.write_subcellular_location_tidy(
    "out/subcellular_location",
    should_write_manifest=True,
)
```

Subcellular location tidy output writes `data.parquet`. It preserves the
UniProtKB `CC   -!- SUBCELLULAR LOCATION:` comment text and ECO evidence; it
does not infer GO cellular-component terms.

## Development

- Project architecture, test standards, benchmarks, and archived plans are
  indexed in [docs/README.md](docs/README.md).
- `pdm run format`
- `pdm run lint`
- `pdm run typecheck`
- `pdm run test`
- `pdm run precommit` runs the complete local gate in that order.
- `PYTHONPATH=src python scripts/benchmark_stringdb.py`

## Release

- GitHub Actions now provides:
  - `.github/workflows/py-ci.yml` for test-and-build checks on push and pull request
  - `.github/workflows/publish.yml` for tag-triggered PyPI publishing
- Release tags must be canonical PEP 440 versions such as `0.1.1`
- The publish workflow expects PyPI trusted publishing to be configured for the `pypi` environment
