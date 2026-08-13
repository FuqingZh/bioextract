from __future__ import annotations

import os
import re
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.errors import CapabilityError, IntegrityError
from bioextract.wikipathways import WikiPathwaysDatabase


def write_gmt(path: Path, *rows: str) -> Path:
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_wikipathways_fixture(tmp_path: Path) -> Path:
    file_gmt = tmp_path / "wikipathways.gmt"
    file_gmt.write_text(
        "\n".join(
            [
                "Glutathione metabolism%WikiPathways_20260510%WP100%Homo sapiens\thttps://www.wikipathways.org/instance/WP100\t2687\t2678\t2678",
                "Alanine and aspartate metabolism%WikiPathways_20260510%WP106%Homo sapiens\thttps://www.wikipathways.org/instance/WP106\t2806\t435",
                "Mouse pathway%WikiPathways_20260510%WP1%Mus musculus\thttps://www.wikipathways.org/instance/WP1\t123",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return file_gmt


def test_from_gmt_resolves_literal_sequence_and_glob_deterministically(
    tmp_path: Path,
) -> None:
    file_b = write_gmt(
        tmp_path / "b.gmt",
        "Mouse pathway%WikiPathways_20260510%WP2%Mus musculus\thttps://example/WP2\t2",
    )
    file_a = write_gmt(
        tmp_path / "a.gmt",
        "Human pathway%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )

    literal = WikiPathwaysDatabase.from_gmt(file_a, glob=False)
    assert literal.snapshot.files_gmt == (file_a.resolve(),)
    sequence = WikiPathwaysDatabase.from_gmt([file_b, file_a], glob=False)
    globbed = WikiPathwaysDatabase.from_gmt(str(tmp_path / "*.gmt"))
    assert (
        sequence.snapshot.files_gmt
        == globbed.snapshot.files_gmt
        == (
            file_a.resolve(),
            file_b.resolve(),
        )
    )
    assert sequence.pathways().collect()["wiki_pathways_id"].to_list() == ["WP1", "WP2"]


def test_from_gmt_glob_false_treats_patterns_literally(tmp_path: Path) -> None:
    write_wikipathways_fixture(tmp_path)

    with pytest.raises(FileNotFoundError, match="file not found"):
        WikiPathwaysDatabase.from_gmt(tmp_path / "*.gmt", glob=False)


def test_from_gmt_double_star_glob_is_recursive(tmp_path: Path) -> None:
    nested = tmp_path / "nested" / "snapshot"
    nested.mkdir(parents=True)
    file_gmt = write_gmt(
        nested / "pathways.gmt",
        "Nested%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )

    db = WikiPathwaysDatabase.from_gmt(str(tmp_path / "**" / "*.gmt"))
    assert db.snapshot.files_gmt == (file_gmt.resolve(),)


@pytest.mark.parametrize("source", [[], ()])
def test_from_gmt_rejects_empty_source(
    source: list[Path] | tuple[Path, ...],
) -> None:
    with pytest.raises(ValueError, match="at least one path"):
        WikiPathwaysDatabase.from_gmt(source)


def test_from_gmt_rejects_empty_scalar_path() -> None:
    with pytest.raises(ValueError, match="paths must be non-empty"):
        WikiPathwaysDatabase.from_gmt("")


def test_from_gmt_rejects_unmatched_missing_directory_and_non_file(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "directory"
    directory.mkdir()
    fifo = tmp_path / "pipe.gmt"
    os.mkfifo(fifo)

    with pytest.raises(FileNotFoundError, match="matched no files"):
        WikiPathwaysDatabase.from_gmt(tmp_path / "*.missing")
    with pytest.raises(FileNotFoundError, match="file not found"):
        WikiPathwaysDatabase.from_gmt(tmp_path / "missing.gmt", glob=False)
    with pytest.raises(ValueError, match="not a file"):
        WikiPathwaysDatabase.from_gmt(directory, glob=False)
    with pytest.raises(ValueError, match="not a file"):
        WikiPathwaysDatabase.from_gmt(fifo, glob=False)


def test_from_gmt_rejects_duplicate_physical_files(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    alias = tmp_path / "alias.gmt"
    alias.symlink_to(file_gmt)
    hard_link = tmp_path / "hard-link.gmt"
    hard_link.hardlink_to(file_gmt)

    duplicate_sources = (
        [file_gmt, file_gmt],
        [str(tmp_path / "*.gmt"), file_gmt],
        [file_gmt, alias],
        [file_gmt, hard_link],
    )
    for source in duplicate_sources:
        with pytest.raises(ValueError, match="duplicate physical file"):
            WikiPathwaysDatabase.from_gmt(source)


def test_from_gmt_requires_common_collection(
    tmp_path: Path,
) -> None:
    file_a = write_gmt(
        tmp_path / "a.gmt",
        "A%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )
    file_b = write_gmt(
        tmp_path / "b.gmt",
        "B%WikiPathways_20260511%WP2%Mus musculus\thttps://example/WP2\t2",
    )

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="common Collection"
    ):
        WikiPathwaysDatabase.from_gmt([file_a, file_b]).pathways().collect()


def test_from_gmt_rejects_empty_file_alongside_valid_file(tmp_path: Path) -> None:
    valid = write_gmt(
        tmp_path / "valid.gmt",
        "A%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )
    empty = tmp_path / "empty.gmt"
    empty.write_text("\n \t\n", encoding="utf-8")

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError),
        match=rf"at least one non-empty pathway record: path={re.escape(str(empty))}",
    ):
        WikiPathwaysDatabase.from_gmt([valid, empty]).pathways().collect()


def test_from_gmt_rejects_single_empty_file(tmp_path: Path) -> None:
    empty = tmp_path / "empty.gmt"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError),
        match=rf"at least one non-empty pathway record: path={re.escape(str(empty))}",
    ):
        WikiPathwaysDatabase.from_gmt(empty).pathways().collect()


@pytest.mark.parametrize("collection", ["Reactome_20260510", "WikiPathways_"])
def test_from_gmt_rejects_invalid_collection_version(
    tmp_path: Path,
    collection: str,
) -> None:
    file_gmt = write_gmt(
        tmp_path / "invalid-collection.gmt",
        f"A%{collection}%WP1%Homo sapiens\thttps://example/WP1\t1",
    )

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="WikiPathways GMT Collection"
    ):
        WikiPathwaysDatabase.from_gmt(file_gmt).pathways().collect()


def test_from_gmt_extracts_one_common_version_from_collection(
    tmp_path: Path,
) -> None:
    file_a = write_gmt(
        tmp_path / "a.gmt",
        "A%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )
    file_b = write_gmt(
        tmp_path / "b.gmt",
        "B%WikiPathways_20260510%WP2%Mus musculus\thttps://example/WP2\t2",
    )

    pathways = WikiPathwaysDatabase.from_gmt([file_a, file_b]).pathways().collect()
    assert pathways["collection"].unique().to_list() == ["WikiPathways_20260510"]
    assert pathways["version"].unique().to_list() == ["20260510"]


@pytest.mark.parametrize("split_files", [False, True])
def test_from_gmt_rejects_duplicate_ids_within_or_across_files(
    tmp_path: Path,
    split_files: bool,
) -> None:
    row_a = "A%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/A\t1"
    row_b = "B%WikiPathways_20260510%WP1%Mus musculus\thttps://example/B\t2"
    sources = (
        [
            write_gmt(tmp_path / "a.gmt", row_a),
            write_gmt(tmp_path / "b.gmt", row_b),
        ]
        if split_files
        else write_gmt(tmp_path / "both.gmt", row_a, row_b)
    )

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError),
        match="wiki_pathways_id must be unique",
    ):
        WikiPathwaysDatabase.from_gmt(sources).pathways().collect()


def test_extract_pathway_term_frames_and_species_filter(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt).with_species("Homo sapiens")

    assert db.pathways().collect().to_dicts() == [
        {
            "wiki_pathways_id": "WP100",
            "pathway_name": "Glutathione metabolism",
            "species": "Homo sapiens",
            "collection": "WikiPathways_20260510",
            "version": "20260510",
            "url": "https://www.wikipathways.org/instance/WP100",
            "gene_count": 2,
        },
        {
            "wiki_pathways_id": "WP106",
            "pathway_name": "Alanine and aspartate metabolism",
            "species": "Homo sapiens",
            "collection": "WikiPathways_20260510",
            "version": "20260510",
            "url": "https://www.wikipathways.org/instance/WP106",
            "gene_count": 2,
        },
    ]
    assert db.pathway_genes().collect().to_dicts() == [
        {"wiki_pathways_id": "WP100", "gene_id": "2678"},
        {"wiki_pathways_id": "WP100", "gene_id": "2687"},
        {"wiki_pathways_id": "WP106", "gene_id": "2806"},
        {"wiki_pathways_id": "WP106", "gene_id": "435"},
    ]
    assert db.pathway_names().collect().columns == [
        "wiki_pathways_id",
        "pathway_name",
        "species",
        "collection",
        "version",
        "url",
    ]


def test_single_and_grouped_selection(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt).with_species("Homo sapiens")

    selection = db.select_ids(["2687", " 435 ", "MISSING", ""])
    assert selection.mappings().collect().to_dicts() == [
        {
            "input_id": "2687",
            "gene_id": "2687",
            "wiki_pathways_id": "WP100",
            "pathway_name": "Glutathione metabolism",
            "species": "Homo sapiens",
            "url": "https://www.wikipathways.org/instance/WP100",
        },
        {
            "input_id": "435",
            "gene_id": "435",
            "wiki_pathways_id": "WP106",
            "pathway_name": "Alanine and aspartate metabolism",
            "species": "Homo sapiens",
            "url": "https://www.wikipathways.org/instance/WP106",
        },
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [{"input_id": "MISSING"}]

    grouped = db.select_groups({"A": ["2687"], "B": ["2687", "MISSING"]})
    df_grouped_mapping = grouped.mappings().collect()
    assert grouped.mappings().collect().equals(df_grouped_mapping)
    assert df_grouped_mapping.columns[0] == "group_id"
    assert df_grouped_mapping.height == 2
    df_unmatched = grouped.unmatched_ids().collect()
    assert grouped.unmatched_ids().collect().equals(df_unmatched)
    assert df_unmatched.to_dicts() == [{"group_id": "B", "input_id": "MISSING"}]


def test_build_tidy_writes_duckdb_without_sidecar(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt).with_species("Homo sapiens")

    tidy = db.build_tidy()
    assert set(tidy.frames) == {"pathway", "term2gene", "term2name"}
    result = db.write_duckdb(tmp_path / "wikipathways.duckdb")
    assert result.tables == ("pathway", "pathway_gene")
    assert not (tmp_path / "manifest.json").exists()


def test_species_filter_keeps_all_file_provenance(tmp_path: Path) -> None:
    file_human = write_gmt(
        tmp_path / "human.gmt",
        "Human%WikiPathways_20260510%WP1%Homo sapiens\thttps://example/WP1\t1",
    )
    file_mouse = write_gmt(
        tmp_path / "mouse.gmt",
        "Mouse%WikiPathways_20260510%WP2%Mus musculus\thttps://example/WP2\t2",
    )
    db = WikiPathwaysDatabase.from_gmt([file_mouse, file_human]).with_species(
        "Homo sapiens"
    )

    assert db.pathways().collect()["wiki_pathways_id"].to_list() == ["WP1"]
    path_out = tmp_path / "wikipathways.duckdb"
    db.write_duckdb(path_out)
    with duckdb.connect(str(path_out), read_only=True) as connection:
        sources = connection.execute(
            "SELECT logical_name, display_path "
            "FROM _bioextract.source_file ORDER BY logical_name"
        ).fetchall()
        release_metadata = dict(
            connection.execute(
                "SELECT key, value FROM _bioextract.metadata "
                "WHERE key IN ('bioextract.release_version', "
                "'bioextract.release_version_source')"
            ).fetchall()
        )
    assert sources == [
        ("pathway_gmt_001", str(file_human.resolve())),
        ("pathway_gmt_002", str(file_mouse.resolve())),
    ]
    assert release_metadata == {
        "bioextract.release_version": "20260510",
        "bioextract.release_version_source": "official_metadata",
    }


def test_duckdb_reopen_matches_source_extraction_selection_and_native_sql(
    tmp_path: Path,
) -> None:
    source = WikiPathwaysDatabase.from_gmt(
        write_wikipathways_fixture(tmp_path)
    ).with_species("Homo sapiens")
    publication = tmp_path / "wikipathways.duckdb"
    source.write_duckdb(publication)
    reopened = WikiPathwaysDatabase.from_duckdb(publication)

    assert reopened.pathways().collect().equals(source.pathways().collect())
    assert reopened.pathway_genes().collect().equals(source.pathway_genes().collect())
    assert reopened.pathway_names().collect().equals(source.pathway_names().collect())
    assert (
        reopened.select_ids(["2687", "MISSING"])
        .mappings()
        .collect()
        .equals(source.select_ids(["2687", "MISSING"]).mappings().collect())
    )
    assert (
        reopened.select_groups({"case": ["2687"], "control": ["435", "MISSING"]})
        .mappings()
        .collect()
        .equals(
            source.select_groups({"case": ["2687"], "control": ["435", "MISSING"]})
            .mappings()
            .collect()
        )
    )

    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        assert first.sql("SELECT count(*) FROM pathway").fetchone() == (2,)
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_with_species_scopes_pathway_genes_before_entrez_matching(
    tmp_path: Path,
) -> None:
    source_file = write_gmt(
        tmp_path / "same-gene.gmt",
        "Human pathway%WikiPathways_20260510%WPH%Homo sapiens	https://example/WPH	7",
        "Mouse pathway%WikiPathways_20260510%WPM%Mus musculus	https://example/WPM	7",
    )
    database = WikiPathwaysDatabase.from_gmt(source_file)
    human = database.with_species(" Homo sapiens ")

    assert database.species is None
    assert human.species == "Homo sapiens"
    assert isinstance(human.pathways(), pl.LazyFrame)
    assert human.pathways().collect()["wiki_pathways_id"].to_list() == ["WPH"]
    assert human.pathway_genes().collect().to_dicts() == [
        {"wiki_pathways_id": "WPH", "gene_id": "7"}
    ]
    assert human.select_ids(["7"]).mappings().collect().select(
        "input_id", "wiki_pathways_id"
    ).to_dicts() == [{"input_id": "7", "wiki_pathways_id": "WPH"}]
    assert human.select_ids(["8"]).unmatched_ids().collect().to_dicts() == [
        {"input_id": "8"}
    ]

    scoped_publication = tmp_path / "same-gene-human.duckdb"
    human.write_duckdb(scoped_publication)
    with duckdb.connect(str(scoped_publication), read_only=True) as connection:
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata WHERE key = 'bioextract.scope'"
        ).fetchone() == ('{"species":"Homo sapiens"}',)
    reopened_scoped = WikiPathwaysDatabase.from_duckdb(scoped_publication)
    assert reopened_scoped.species == "Homo sapiens"
    with pytest.raises(CapabilityError, match="scoped to"):
        reopened_scoped.with_species("Mus musculus")

    publication = tmp_path / "same-gene.duckdb"
    database.write_duckdb(publication)
    reopened_human = WikiPathwaysDatabase.from_duckdb(publication).with_species(
        "Homo sapiens"
    )
    assert reopened_human.pathways().collect()["wiki_pathways_id"].to_list() == ["WPH"]
    assert reopened_human.select_ids(["7"]).mappings().collect().select(
        "input_id", "wiki_pathways_id"
    ).to_dicts() == [{"input_id": "7", "wiki_pathways_id": "WPH"}]


def test_duckdb_reopen_validates_bounded_publication_contract(tmp_path: Path) -> None:
    source_file = write_wikipathways_fixture(tmp_path)

    def publish(name: str) -> Path:
        path = tmp_path / name
        WikiPathwaysDatabase.from_gmt(source_file).write_duckdb(path)
        return path

    wrong_profile = publish("wrong-profile.duckdb")
    with duckdb.connect(str(wrong_profile)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='forged-v1' "
            "WHERE key='bioextract.source_schema_profile'"
        )
    with pytest.raises(IntegrityError, match="source schema profile"):
        WikiPathwaysDatabase.from_duckdb(wrong_profile)

    inventory_drift = publish("inventory-drift.duckdb")
    with duckdb.connect(str(inventory_drift)) as connection:
        connection.execute("CREATE VIEW unrecorded AS SELECT * FROM pathway")
    with pytest.raises(IntegrityError, match="table/view inventory"):
        WikiPathwaysDatabase.from_duckdb(inventory_drift)

    schema_drift = publish("schema-drift.duckdb")
    with duckdb.connect(str(schema_drift)) as connection:
        connection.execute(
            "ALTER TABLE pathway ALTER gene_count TYPE VARCHAR USING gene_count::VARCHAR"
        )
    with pytest.raises(IntegrityError, match="table schema"):
        WikiPathwaysDatabase.from_duckdb(schema_drift)


def test_reopened_capabilities_and_cached_terminals_recheck_identity(
    tmp_path: Path,
) -> None:
    source_file = write_wikipathways_fixture(tmp_path)
    publication = tmp_path / "wikipathways.duckdb"
    source = WikiPathwaysDatabase.from_gmt(source_file)
    source.write_duckdb(publication)
    reopened = WikiPathwaysDatabase.from_duckdb(publication)
    selection = reopened.select_ids(["2687"])
    selection.mappings().collect()

    with pytest.raises(CapabilityError, match="GMT source handle"):
        reopened.write_duckdb(tmp_path / "copy.duckdb")
    with pytest.raises(CapabilityError, match="GMT source handle"):
        reopened.build_tidy().write_duckdb(tmp_path / "copy.duckdb")
    with pytest.raises(CapabilityError, match="from_duckdb"):
        source.connect()

    replacement_dir = tmp_path / "replacement"
    replacement_dir.mkdir()
    replacement = replacement_dir / publication.name
    WikiPathwaysDatabase.from_gmt(source_file).write_duckdb(replacement)
    replacement.replace(publication)

    with pytest.raises((IntegrityError, pl.exceptions.ComputeError), match="replaced"):
        selection.mappings().collect()
    with pytest.raises((IntegrityError, pl.exceptions.ComputeError), match="replaced"):
        reopened.pathways().collect()
    with pytest.raises(IntegrityError, match="replaced"):
        reopened.connect()


def test_from_gmt_rejects_missing_and_malformed_files(
    tmp_path: Path,
) -> None:
    write_wikipathways_fixture(tmp_path)

    with pytest.raises(FileNotFoundError):
        WikiPathwaysDatabase.from_gmt(tmp_path / "missing.gmt")

    file_bad = tmp_path / "bad.gmt"
    file_bad.write_text("bad\thttps://example.org\t1\n", encoding="utf-8")
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError),
        match="four '%' separated fields",
    ):
        WikiPathwaysDatabase.from_gmt(file_bad).pathways().collect()
