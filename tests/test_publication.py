from __future__ import annotations

import inspect
from pathlib import Path

import duckdb
import polars as pl
import pytest

import bioextract.chebi as chebi
import bioextract.eggnog as eggnog
import bioextract.go as go
import bioextract.interpro as interpro
import bioextract.kegg as kegg
import bioextract.omnipath as omnipath
import bioextract.reactome as reactome
import bioextract.rhea as rhea
import bioextract.stringdb as stringdb
import bioextract.uniprot as uniprot
import bioextract.wikipathways as wikipathways
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.stringdb.stringdb import StringSelection


def _dataset(tmp_path: Path, *, relation_count: int = 1) -> TidyDataset:
    source = tmp_path / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    frames = {
        f"relation_{index}": pl.DataFrame({"id": [f"T{index}"]}).lazy()
        for index in range(relation_count)
    }
    return TidyDataset(
        frames=frames,
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=tuple(
            TidyAsset(
                path=f"relation_{index}.parquet",
                kind="canonical",
                frame_name=f"relation_{index}",
            )
            for index in range(relation_count)
        ),
        resource_name="example",
        release_version="2026-07-29",
    )


def test_parquet_publication_embeds_provenance_without_sidecar(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.parquet"
    result = _dataset(tmp_path).write_parquet(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    metadata = dict(
        duckdb.connect()
        .execute(
            "SELECT key, value FROM parquet_kv_metadata(?) "
            "WHERE CAST(key AS VARCHAR) LIKE 'bioextract.%'",
            [str(path)],
        )
        .fetchall()
    )
    assert metadata[b"bioextract.resource_name"] == b"example"
    assert metadata[b"bioextract.resource_schema_version"] == b"example-v1"
    assert metadata[b"bioextract.release_version"] == b"2026-07-29"


@pytest.mark.parametrize(
    ("release_version", "release_version_source", "message"),
    [
        (" ", None, "release_version must be non-empty"),
        ("2026_01", "filename", "caller or official_metadata"),
        (None, "caller", "requires release_version"),
    ],
)
def test_publication_rejects_invalid_release_provenance(
    tmp_path: Path,
    release_version: str | None,
    release_version_source: str | None,
    message: str,
) -> None:
    dataset = _dataset(tmp_path)
    dataset.release_version = release_version
    dataset.release_version_source = release_version_source
    with pytest.raises(ValueError, match=message):
        dataset.write_parquet(tmp_path / "invalid-release.parquet")


def test_duckdb_publication_has_internal_provenance_schema(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.duckdb"
    result = _dataset(tmp_path, relation_count=2).write_duckdb(path)

    assert result.tables == ("relation_0", "relation_1")
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata_tables = {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = '_bioextract'"
            ).fetchall()
        }
        assert metadata_tables == {
            "column_mapping",
            "metadata",
            "source_file",
            "table_info",
            "validation_issue",
        }
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key = 'bioextract.metadata_schema_version'"
        ).fetchone() == ("3",)
        assert connection.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone() == (0,)
        resource_row = connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key = 'bioextract.resource_name'"
        ).fetchone()
        assert resource_row is not None
        assert resource_row[0] == "example"
        assert connection.execute(
            "SELECT table_name, table_role, row_count "
            "FROM _bioextract.table_info ORDER BY table_name"
        ).fetchall() == [
            ("relation_0", "canonical", 1),
            ("relation_1", "canonical", 1),
        ]
        mapping_row = connection.execute(
            "SELECT count(*) FROM _bioextract.column_mapping"
        ).fetchone()
        assert mapping_row is not None
        assert mapping_row[0] == 0


def test_canonical_publication_normalizes_derived_columns_and_records_mapping(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("UniProtId\nP12345\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={
            "protein_pathway": pl.DataFrame(
                {"UniProtId": ["P12345"], "ReactomePathwayId": ["R-HSA-1"]}
            ).lazy()
        },
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=(
            TidyAsset(
                "protein_pathway.parquet",
                "canonical",
                "protein_pathway",
            ),
        ),
    )

    file_parquet = tmp_path / "example.parquet"
    dataset.write_parquet(file_parquet)
    assert pl.read_parquet(file_parquet).columns == [
        "uniprot_id",
        "reactome_pathway_id",
    ]

    file_duckdb = tmp_path / "example.duckdb"
    dataset.write_duckdb(file_duckdb)
    with duckdb.connect(str(file_duckdb), read_only=True) as connection:
        assert [
            row[1]
            for row in connection.execute(
                "PRAGMA table_info('protein_pathway')"
            ).fetchall()
        ] == ["uniprot_id", "reactome_pathway_id"]
        assert connection.execute(
            "SELECT source_column, output_column "
            "FROM _bioextract.column_mapping ORDER BY source_column"
        ).fetchall() == [
            ("ReactomePathwayId", "reactome_pathway_id"),
            ("UniProtId", "uniprot_id"),
        ]


def test_official_headers_receive_only_required_duckdb_mapping(
    tmp_path: Path,
) -> None:
    source = tmp_path / "official.tsv"
    source.write_text("Name\tname\nA\tB\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={"official": pl.DataFrame({"Name": ["A"], "name": ["B"]}).lazy()},
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="official-v1",
        source_schema_profile="official-source-v1",
        build_id_prefix="official",
        assets=(TidyAsset("official.parquet", "canonical", "official"),),
    )

    path = tmp_path / "official.duckdb"
    dataset.write_duckdb(
        path,
        preserve_source_headers={"official"},
    )

    with duckdb.connect(str(path), read_only=True) as connection:
        assert [
            row[1]
            for row in connection.execute("PRAGMA table_info('official')").fetchall()
        ] == ["Name", "name_2"]
        assert connection.execute(
            "SELECT source_column, output_column, reason "
            "FROM _bioextract.column_mapping"
        ).fetchall() == [("name", "name_2", "case_insensitive_collision")]


def test_publication_is_atomic_and_requires_snake_case_tables(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.parquet"
    _dataset(tmp_path).write_parquet(path)
    original = path.read_bytes()

    with pytest.raises(FileExistsError):
        _dataset(tmp_path).write_parquet(path)
    assert path.read_bytes() == original

    dataset = _dataset(tmp_path, relation_count=2)
    with pytest.raises(ValueError, match="snake_case"):
        dataset.write_duckdb(
            tmp_path / "bad.duckdb",
            table_names={"relation_0": "relation-zero"},
        )


@pytest.mark.parametrize("container", ["parquet", "duckdb"])
def test_failed_replacement_preserves_existing_publication(
    tmp_path: Path,
    container: str,
) -> None:
    path = tmp_path / f"example.{container}"
    if container == "parquet":
        _dataset(tmp_path).write_parquet(path)
    else:
        _dataset(tmp_path).write_duckdb(path)
    original = path.read_bytes()

    source = tmp_path / "bad-source.tsv"
    source.write_text("value\nbad\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={
            "relation": pl.DataFrame({"value": ["bad"]})
            .lazy()
            .select(pl.col("value").cast(pl.Int64))
        },
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="bad-v1",
        source_schema_profile="bad-source-v1",
        build_id_prefix="bad",
        assets=(TidyAsset("relation.parquet", "canonical", "relation"),),
    )
    with pytest.raises(pl.exceptions.InvalidOperationError):
        if container == "parquet":
            dataset.write_parquet(path, if_exists="replace")
        else:
            dataset.write_duckdb(path, if_exists="replace")

    assert path.read_bytes() == original


@pytest.mark.parametrize(
    ("module", "old_name"),
    [
        (go, "GoDb"),
        (kegg, "KeggDb"),
        (reactome, "ReactomeDb"),
        (wikipathways, "WikiPathwaysDb"),
        (eggnog, "EggnogDb"),
        (interpro, "InterProDb"),
        (uniprot, "UniprotDb"),
        (stringdb, "StringDb"),
        (omnipath, "OmniPathDb"),
        (rhea, "RheaDb"),
    ],
)
def test_legacy_database_type_aliases_are_not_exported(
    module: object,
    old_name: str,
) -> None:
    assert not hasattr(module, old_name)


def test_removed_legacy_writer_apis_are_not_exported() -> None:
    assert not hasattr(TidyDataset, "write")
    assert not hasattr(StringSelection, "with_score_min")


def test_resource_factories_do_not_expose_limits() -> None:
    factories = (
        chebi.ChEBIDatabase.from_release,
        go.GODatabase.from_obo,
        kegg.KEGGDatabase.from_brite_json,
        reactome.ReactomeDatabase.from_files,
        wikipathways.WikiPathwaysDatabase.from_gmt,
        eggnog.EggNOGDatabase.from_files,
        interpro.InterProDatabase.from_mapping_files,
        uniprot.UniProtDatabase.from_idmapping,
        uniprot.UniProtDatabase.from_knowledgebase,
        stringdb.STRINGDatabase.from_files,
        omnipath.OmniPathDatabase.from_files,
        rhea.RheaDatabase.from_release,
    )
    assert all(
        "limits" not in inspect.signature(factory).parameters for factory in factories
    )


def test_resource_factory_parameter_names_follow_domain_roles() -> None:
    expected = {
        chebi.ChEBIDatabase.from_release: ("source", "chemont_obo"),
        chebi.ChEBIDatabase.from_table_files: (
            "compounds",
            "names",
            "relations",
            "secondary_ids",
            "database_accessions",
            "structures",
            "chemical_data",
            "chemont_obo",
        ),
        chebi.ChEBIDatabase.from_obo: ("path", "sdf", "chemont_obo"),
        chebi.ChEBIDatabase.from_duckdb: ("path",),
        go.GODatabase.from_obo: ("path",),
        kegg.KEGGDatabase.from_brite_json: ("path",),
        kegg.KEGGDatabase.from_mapping_files: (
            "uniprot_conversion",
            "gene_ko",
            "gene_pathway",
            "organism_code",
            "gene_list",
            "ncbi_gene_conversion",
        ),
        reactome.ReactomeDatabase.from_files: (
            "uniprot_mapping",
            "pathways",
            "relations",
        ),
        wikipathways.WikiPathwaysDatabase.from_gmt: ("path", "species"),
        eggnog.EggNOGDatabase.from_files: (
            "eggnog_database",
            "cog_functions",
            "temp_dir",
        ),
        interpro.InterProDatabase.from_mapping_files: (
            "protein_to_interpro",
            "interpro_xml",
        ),
        uniprot.UniProtDatabase.from_idmapping: ("path", "release_version"),
        uniprot.UniProtDatabase.from_knowledgebase: (
            "entries",
            "canonical_sequences",
            "isoform_sequences",
            "release_version",
        ),
        uniprot.UniProtDatabase.from_duckdb: ("path",),
        stringdb.STRINGDatabase.from_files: (
            "aliases",
            "links",
            "rank_by_source",
            "release_version",
        ),
        omnipath.OmniPathDatabase.from_files: ("enzsub", "interactions"),
        rhea.RheaDatabase.from_reaction_files: (
            "rdf",
            "directions",
            "relationships",
            "obsolete_reactions",
            "reaction_smiles",
        ),
        rhea.RheaDatabase.from_compound_files: (
            "sdf",
            "chebi_names",
            "chebi_ph7_3_mapping",
        ),
        rhea.RheaDatabase.from_cross_reference_files: (
            "xrefs",
            "uniprot_sprot",
            "uniprot_trembl",
        ),
        rhea.RheaDatabase.from_release: ("source",),
        rhea.RheaDatabase.from_duckdb: ("path",),
    }
    assert {
        factory: tuple(inspect.signature(factory).parameters) for factory in expected
    } == expected
