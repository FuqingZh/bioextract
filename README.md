# bioextract

Stable, provenance-aware domain access to official biological database
snapshots.

`bioextract` hides resource-specific file layouts, identifier rules,
hierarchies, directions, and repeated joins. Callers provide local snapshot
files; the library neither downloads resources nor knows which application
will consume them.

## Architecture

The domain contract is primary. Storage is an execution strategy:

- consume an official indexed or inexpensive representation directly when it
  is already fit for use;
- publish one independent relation as Parquet;
- publish related relations that are normally queried together as one DuckDB;
- keep ordinary filtering, sorting, grouping, and SQL in Polars or DuckDB;
- add convenience methods only when they encode resource-owned ID resolution,
  relationship traversal, grouping, or unmatched-ID accounting.

Canonical writers always require an explicit `path`. They validate into a
staging file and atomically publish only after success. Parquet provenance is
embedded in footer metadata; DuckDB provenance lives in the `_bioextract`
schema. No sidecar manifest is required.

Read the
[Domain Access Architecture](docs/architecture/20260729-v1.0-domain-access-architecture.md)
before adding a resource, public query method, or output format.

## Install

```console
pip install bioextract
```

## ChEBI and ChemOnt

FULL OBO supplies the canonical compound, identifier, name, cross-reference,
property, and relation schema. SDF only supplements molfile records; optional
ChemOnt remains a separate `chemont_*` graph in the same container:

```python
from bioextract.chebi import ChEBIDatabase

result = ChEBIDatabase.from_release(
    "chebi/database/2026-07-07/raw",
    chemont_obo="ChemOnt_2_1.obo.zip",
).write_duckdb("out/chebi.duckdb")

print(result.tables)
```

Open the publication for stable domain extraction or unrestricted native
read-only SQL:

```python
database = ChEBIDatabase.from_duckdb("out/chebi.duckdb")
selection = database.select_compounds(
    ["CHEBI:15377", "CHEBI:10743"],
    namespace="chebi",
)

df_compounds = selection.extract_compounds()
df_names = selection.extract_names()
df_relations = selection.extract_relations()
df_unmatched = selection.extract_unmatched_ids()

with database.connect() as connection:
    prefixes = connection.execute(
        "SELECT DISTINCT source_prefix FROM compound_cross_reference"
    ).fetchall()
```

External cross-references use the official prefix directly as `namespace`,
such as `kegg.compound` or `hmdb`. Public shared IDs are complete
`CHEBI:<number>` CURIEs. Use explicit TSV files only for partial source builds;
plain, gzip, zip, and tar inputs are detected internally where applicable.

## Rhea

Build one query-ready database from a complete extracted release or archive:

```python
from bioextract.rhea import RheaDatabase

result = RheaDatabase.from_release("rhea-release.zip").write_duckdb(
    "out/rhea.duckdb"
)
print(result.tables)
```

Focused constructors accept incomplete capabilities while retaining the same
DuckDB container:

```python
from bioextract.rhea import RheaDatabase

RheaDatabase.from_reaction_files(
    rdf="rhea.rdf.gz",
    directions="rhea-directions.tsv",
    relationships="rhea-relationships.tsv",
).write_duckdb("out/rhea.duckdb")

RheaDatabase.from_cross_reference_files(
    xrefs="rhea2xrefs.tsv",
    uniprot_sprot="rhea2uniprot_sprot.tsv",
    uniprot_trembl="rhea2uniprot_trembl.tsv.gz",
).write_duckdb("out/rhea_xrefs.duckdb")
```

Open a published database and select reactions through any one supported
official namespace:

```python
database = RheaDatabase.from_duckdb("out/rhea.duckdb")
selection = database.select_reactions(
    ["CHEBI:15377", "CHEBI:16474"],
    namespace="chebi",
)

df_matches = selection.extract_matches()
df_reactions = selection.extract_reactions()
df_participants = selection.extract_participants()
df_cross_references = selection.extract_cross_references()
df_unmatched = selection.extract_unmatched_ids()
```

`select_reactions()` and `select_groups()` are deferred domain query plans;
their `extract_*()` terminals return eager Polars `DataFrame` objects.
Participant output retains the exact Rhea ID, master ID, direction, side, and
compound fields. ChEBI fields are complete `CHEBI:<number>` CURIEs and can be
equality-joined to a ChEBI publication without prefix construction or casts.
`DirectionalRole` is populated only for `LR` and `RL`;
undefined and bidirectional reactions retain null rather than inventing a
substrate/product orientation.

See the [Rhea architecture](docs/architecture/rhea-db.md) for direction,
hierarchy, table, and provenance contracts.

## GO

GO is a multi-relation ontology and is published as one DuckDB:

```python
from bioextract.go import GODatabase

go = GODatabase.from_obo("go-basic.obo")
df_terms = go.select_terms(subset_id="goslim_generic")
df_subcell = go.extract_subcell()
result = go.write_duckdb("out/go.duckdb")
```

Tables include `term`, `term_relation`, `term_synonym`, `term_xref`,
`term_alternate_id`, `term_ancestor`, and `term_depth`.

## KEGG

An independent KEGG mapping or BRITE relation is published as Parquet:

```python
from bioextract.kegg import KEGGDatabase

KEGGDatabase.from_brite_json("br08901.json").write_parquet(
    "out/kegg.parquet"
)
```

When multiple KEGG products share a directory, use the smallest useful
qualifier, such as `kegg_gene_annotation.parquet`.

A compound/reaction/enzyme/module snapshot is a multi-relation metabolic
publication:

```python
database = KEGGDatabase.from_metabolic_release("kegg/metabolic/2026-07")
database.write_duckdb("out/kegg.duckdb")

published = KEGGDatabase.from_duckdb("out/kegg.duckdb")
selection = published.select_ids(["CHEBI:15377"], namespace="chebi")

df_reactions = selection.extract_reactions()
df_pathways = selection.extract_pathway_memberships()
df_unmatched = selection.extract_unmatched_ids()

with published.connect() as connection:
    relation = connection.sql(
        """
        SELECT reaction_id, count(*) AS participant_count
        FROM reaction_participant
        GROUP BY reaction_id
        """
    )
```

The domain API supplies reaction-centered traversal and input lineage.
`connect()` exposes the same validated publication as caller-owned, native
read-only DuckDB SQL.

## Reactome and WikiPathways

Pathway entities and membership relations are published together:

```python
from bioextract.reactome import ReactomeDatabase
from bioextract.wikipathways import WikiPathwaysDatabase

ReactomeDatabase.from_files(
    uniprot_mapping="UniProt2Reactome.txt",
    pathways="ReactomePathways.txt",
    relations="ReactomePathwaysRelation.txt",
).write_duckdb("out/reactome.duckdb")

WikiPathwaysDatabase.from_gmt(
    "wikipathways-20260510-gmt-*.gmt",
    species="Homo sapiens",
).write_duckdb("out/wikipathways.duckdb")
```

Selection methods such as `select_ids()` retain unmatched inputs and hide
resource-specific mapping joins. They do not calculate enrichment statistics.
WikiPathways glob expansion is enabled by default; pass `glob=False` for a
literal path or sequence. The constructor validates one Collection and Version
and unique pathway IDs across the complete resolved file set before applying
the optional species row filter.

## eggNOG

The canonical protein-to-orthologous-group mapping is one Parquet:

```python
from bioextract.eggnog import EggNOGDatabase

EggNOGDatabase.from_sqlite(
    "eggnog.db.gz",
    cog_functions="cog-24.fun.tab",
).write_parquet("out/eggnog.parquet")
```

The official SQLite representation is consumed directly during extraction;
publication is optional and exists for repeated analytical scans.

## InterPro and Pfam

The independent InterPro mapping is Parquet. The related Pfam term, xref, and
protein-term relations share one DuckDB:

```python
from bioextract.interpro import InterProDatabase

db = InterProDatabase.from_mapping_files(
    protein_to_interpro="108.0/raw/protein2ipr.dat.gz",
    interpro_xml="108.0/raw/interpro.xml.gz",
)
db.write_parquet("out/interpro.parquet")
db.write_duckdb("out/interpro_pfam.duckdb", config="pfam")
```

## UniProt

UniProt idmapping remains a separate lazy Parquet product:

```python
from bioextract.uniprot import UniProtDatabase

UniProtDatabase.from_idmapping(
    "idmapping_selected.tab.gz",
    release_version="2026_01",
).write_parquet(
    "out/uniprot_idmapping.parquet",
    taxon_ids=["9606", "10090"],
)
```

Reviewed UniProtKB is a multi-relation DuckDB publication:

```python
UniProtDatabase.from_knowledgebase(
    entries="uniprot_sprot.dat.gz",
    canonical_sequences="uniprot_sprot.fasta.gz",
    isoform_sequences="uniprot_sprot_varsplic.fasta.gz",
    release_version="2026_01",
).write_duckdb("out/uniprot.duckdb")

db = UniProtDatabase.from_duckdb("out/uniprot.duckdb")
proteins = db.select_ids(
    ["P04637"],
    namespace="uniprot",
    taxon_ids=["9606"],
).extract_proteins()
with db.connect() as connection:
    relation_count = connection.execute(
        "SELECT count(*) FROM protein"
    ).fetchone()[0]
```

Constructor arguments declare source roles, while headers and record grammar
validate their content. Paths never supply release identity. An all-taxid
idmapping export requires `allow_all_taxa=True`.

## STRING

`select_ids()` and `select_groups()` encapsulate alias resolution, unmatched
IDs, group isolation, and edge mapping:

```python
from bioextract.stringdb import STRINGDatabase

selection = (
    STRINGDatabase.from_files(
        aliases="9606.protein.aliases.v12.0.txt.gz",
        links="9606.protein.links.v12.0.txt.gz",
    )
    .select_groups(
        {
            "TumorA": ["TP53", "EGFR"],
            "TumorB": ["CDK2", "TP53"],
        }
    )
    .with_min_combined_score(400)
)

df_mapping = selection.extract_string_mapping()
df_unmapped = selection.extract_unmatched_ids()
df_edges = selection.extract_edges()
```

`combined_score` is a STRING confidence score, not an interaction-strength
measurement.

## OmniPath

```python
from bioextract.omnipath import OmniPathDatabase

selection = (
    OmniPathDatabase.from_files(
        enzsub="enzsub.tsv.gz",
        interactions="interactions.tsv.gz",
    )
    .select_ids(["P31749", "AKT1", "BAD"])
    .with_enzsub()
)

df_enzsub = selection.extract_enzsub()
df_unmapped = selection.extract_unmatched_ids()
```

## Naming and compatibility

Public resource handles use complete `*Database` names, including
`GODatabase`, `ChEBIDatabase`, `RheaDatabase`, `KEGGDatabase`,
`ReactomeDatabase`, `WikiPathwaysDatabase`, `EggNOGDatabase`,
`InterProDatabase`, `UniProtDatabase`, `STRINGDatabase`, and
`OmniPathDatabase`.

There are no abbreviated `*Db` aliases, legacy score-filter names, or
directory writers. Use `with_min_combined_score()`, `write_parquet()`, and
`write_duckdb()` directly.

Table names, view names, and generated columns use singular `snake_case`.
Official two-dimensional source headers are retained unless a minimal
deterministic mapping is required to make them queryable. Any such mapping is
recorded in embedded provenance.

Example names such as `go.duckdb` and `kegg.parquet` are product guidance, not
format standards or compatibility identifiers. Machine identity comes from
embedded metadata.

## Development

- Documentation is indexed in [docs/README.md](docs/README.md).
- `pdm run format`
- `pdm run lint`
- `pdm run typecheck`
- `pdm run test`
- `pdm run precommit` runs the complete local gate in that order.

## Release

- `.github/workflows/py-ci.yml` runs test-and-build checks.
- `.github/workflows/publish.yml` publishes canonical PEP 440 tags.
- PyPI trusted publishing is expected for the `pypi` environment.
