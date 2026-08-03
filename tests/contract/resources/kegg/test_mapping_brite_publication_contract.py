from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import duckdb
import pytest

from bioextract.errors import CapabilityError
from bioextract.kegg import KEGGDatabase


def _mapping_source(tmp_path: Path) -> KEGGDatabase:
    uniprot = tmp_path / "conv_uniprot.tsv"
    gene_ko = tmp_path / "gene_ko.tsv"
    gene_pathway = tmp_path / "gene_pathway.tsv"
    uniprot.write_text("up:P12345\thsa:1\n", encoding="utf-8")
    gene_ko.write_text("hsa:1\tko:K00001\n", encoding="utf-8")
    gene_pathway.write_text("hsa:1\tpath:hsa00010\n", encoding="utf-8")
    return KEGGDatabase.from_mapping_files(
        uniprot_conversion=uniprot,
        gene_ko=gene_ko,
        gene_pathway=gene_pathway,
        organism_code="hsa",
    )


def _brite_source(tmp_path: Path) -> KEGGDatabase:
    source = tmp_path / "brite.json"
    source.write_text(
        json.dumps(
            {
                "name": "hsa00001",
                "children": [
                    {
                        "name": "09100 Metabolism",
                        "children": [
                            {
                                "name": "09101 Carbohydrate metabolism",
                                "children": [
                                    {
                                        "name": "00010 Glycolysis [PATH:hsa00010]",
                                        "children": [
                                            {"name": "hsa:1 GENE\tK00001 enzyme"}
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return KEGGDatabase.from_brite_json(source)


@pytest.mark.parametrize(
    ("source_factory", "scope", "profile", "table"),
    [
        (
            _mapping_source,
            "mapping",
            "kegg-organism-mapping-files-v1",
            "mapping",
        ),
        (_brite_source, "brite", "kegg-brite-json-v1", "pathway"),
    ],
)
def test_single_relation_profiles_publish_metadata_v1_and_reopen(
    tmp_path: Path,
    source_factory: Callable[[Path], KEGGDatabase],
    scope: str,
    profile: str,
    table: str,
) -> None:
    path = tmp_path / f"{scope}.duckdb"
    result = source_factory(tmp_path).write_duckdb(path)

    assert result.tables == (table,)
    reopened = KEGGDatabase.from_duckdb(path)
    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        metadata = dict(
            first.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "1"
        assert metadata["bioextract.resource_name"] == "kegg"
        assert metadata["bioextract.scope"] == scope
        assert metadata["bioextract.source_schema_profile"] == profile
        assert first.execute(
            "SELECT table_name, table_role FROM _bioextract.table_info"
        ).fetchall() == [(table, "canonical")]
        assert first.execute(f"SELECT count(*) FROM {table}").fetchone() == (1,)
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_mapping_publication_preserves_selection_semantics(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    source = _mapping_source(tmp_path)
    expected = source.select_groups(
        {"case": ["P12345", "missing"], "repeat": ["P12345"]},
        namespace="uniprot",
    )
    source.write_duckdb(path)
    actual = KEGGDatabase.from_duckdb(path).select_groups(
        {"case": ["P12345", "missing"], "repeat": ["P12345"]},
        namespace="uniprot",
    )

    assert actual.extract_mapping().equals(expected.extract_mapping())
    assert actual.extract_unmatched_ids().equals(expected.extract_unmatched_ids())


def test_empty_mapping_publication_reopens_with_stable_schema(tmp_path: Path) -> None:
    source = _mapping_source(tmp_path)
    for path in tmp_path.glob("*.tsv"):
        path.write_text("", encoding="utf-8")
    publication = tmp_path / "empty.duckdb"
    source.write_duckdb(publication)

    reopened = KEGGDatabase.from_duckdb(publication).extract_mapping()
    assert reopened.schema == source.extract_mapping().schema
    assert reopened.is_empty()


def test_profile_inventory_and_source_role_mismatches_are_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    _mapping_source(tmp_path).write_duckdb(path)

    with duckdb.connect(str(path)) as connection:
        connection.execute("UPDATE _bioextract.table_info SET table_role='unexpected'")
    with pytest.raises(ValueError, match="table inventory"):
        KEGGDatabase.from_duckdb(path)


def test_atomic_if_exists_preserves_previous_publication(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    source = _mapping_source(tmp_path)
    source.write_duckdb(path)
    before = path.read_bytes()

    with pytest.raises(FileExistsError):
        source.write_duckdb(path)
    assert path.read_bytes() == before


def test_brite_publication_rejects_mapping_selection(tmp_path: Path) -> None:
    path = tmp_path / "brite.duckdb"
    _brite_source(tmp_path).write_duckdb(path)
    database = KEGGDatabase.from_duckdb(path)

    with pytest.raises((CapabilityError, ValueError), match="BRITE"):
        database.select_ids(["P12345"], namespace="uniprot")
