from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

import bioextract.kegg.kegg as kegg_module
from bioextract.kegg import KEGGDatabase


def write_kegg_mapping_fixture(tmp_path: Path) -> dict[str, Path]:
    files = {
        "conv_uniprot": tmp_path / "conv_uniprot.tsv",
        "conv_ncbi_gene": tmp_path / "conv_ncbi_gene.tsv",
        "gene_ko": tmp_path / "gene_ko.tsv",
        "gene_pathway": tmp_path / "gene_pathway.tsv",
        "gene_list": tmp_path / "gene_list.tsv",
    }
    files["conv_uniprot"].write_text(
        "up:P12345\thsa:1\nup:Q9Y243\thsa:2\nup:P12345\thsa:1\n",
        encoding="utf-8",
    )
    files["conv_ncbi_gene"].write_text(
        "ncbi-geneid:101\thsa:1\nncbi-geneid:102\thsa:2\n",
        encoding="utf-8",
    )
    files["gene_ko"].write_text(
        "hsa:1\tko:K00001\nhsa:2\tko:K00002\n",
        encoding="utf-8",
    )
    files["gene_pathway"].write_text(
        "hsa:1\tpath:hsa00010\nhsa:1\tpath:hsa01100\nhsa:2\tpath:hsa04110\n",
        encoding="utf-8",
    )
    files["gene_list"].write_text(
        "hsa:1\tGENE1; alpha description\nhsa:2\tGENE2\n",
        encoding="utf-8",
    )
    return files


def create_mapping_db(files: dict[str, Path]) -> KEGGDatabase:
    return KEGGDatabase.from_mapping_files(
        uniprot_conversion=files["conv_uniprot"],
        ncbi_gene_conversion=files["conv_ncbi_gene"],
        gene_ko=files["gene_ko"],
        gene_pathway=files["gene_pathway"],
        gene_list=files["gene_list"],
        organism_code="hsa",
    )


def test_extract_mapping_normalizes_and_expands_many_to_many(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    df_mapping = db.mappings().collect()

    assert df_mapping.columns == [
        "organism_code",
        "kegg_gene_id",
        "uniprot_id",
        "ncbi_gene_id",
        "ko_id",
        "kegg_pathway_id",
        "pathway_map_id",
        "gene_symbol",
        "gene_description",
    ]
    assert df_mapping.to_dicts() == [
        {
            "organism_code": "hsa",
            "kegg_gene_id": "hsa:1",
            "uniprot_id": "P12345",
            "ncbi_gene_id": "101",
            "ko_id": "K00001",
            "kegg_pathway_id": "hsa00010",
            "pathway_map_id": "map00010",
            "gene_symbol": "GENE1",
            "gene_description": "alpha description",
        },
        {
            "organism_code": "hsa",
            "kegg_gene_id": "hsa:1",
            "uniprot_id": "P12345",
            "ncbi_gene_id": "101",
            "ko_id": "K00001",
            "kegg_pathway_id": "hsa01100",
            "pathway_map_id": "map01100",
            "gene_symbol": "GENE1",
            "gene_description": "alpha description",
        },
        {
            "organism_code": "hsa",
            "kegg_gene_id": "hsa:2",
            "uniprot_id": "Q9Y243",
            "ncbi_gene_id": "102",
            "ko_id": "K00002",
            "kegg_pathway_id": "hsa04110",
            "pathway_map_id": "map04110",
            "gene_symbol": "GENE2",
            "gene_description": None,
        },
    ]


def test_select_ids_supports_input_id_kinds_and_unmapped(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    df_uniprot = (
        db.select_ids(
            ["sp|P12345|GENE1_HUMAN", "MISSING"],
            namespace="uniprot",
        )
        .mappings()
        .collect()
    )
    assert df_uniprot.select(
        "input_id", "input_namespace", "kegg_gene_id"
    ).to_dicts() == [
        {"input_id": "P12345", "input_namespace": "uniprot", "kegg_gene_id": "hsa:1"},
        {"input_id": "P12345", "input_namespace": "uniprot", "kegg_gene_id": "hsa:1"},
    ]
    assert db.select_ids(
        ["P12345", "MISSING"],
        namespace="uniprot",
    ).unmatched_ids().collect().to_dicts() == [{"input_id": "MISSING"}]

    df_ncbi = db.select_ids(["102"], namespace="ncbi_gene").mappings().collect()
    assert df_ncbi.select("input_id", "kegg_gene_id").to_dicts() == [
        {"input_id": "102", "kegg_gene_id": "hsa:2"}
    ]

    df_kegg = db.select_ids(["hsa:1"], namespace="kegg_gene").mappings().collect()
    assert df_kegg.select("input_id", "uniprot_id", "kegg_pathway_id").to_dicts() == [
        {"input_id": "hsa:1", "uniprot_id": "P12345", "kegg_pathway_id": "hsa00010"},
        {"input_id": "hsa:1", "uniprot_id": "P12345", "kegg_pathway_id": "hsa01100"},
    ]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    selection = db.select_groups(
        {"up": ["P12345", "MISSING"], "down": ["Q9Y243"]},
        namespace="uniprot",
    )

    df_mapping = selection.mappings().collect()
    assert df_mapping.columns[:3] == ["group_id", "input_id", "input_namespace"]
    assert df_mapping.select("group_id", "input_id", "kegg_gene_id").to_dicts() == [
        {"group_id": "down", "input_id": "Q9Y243", "kegg_gene_id": "hsa:2"},
        {"group_id": "up", "input_id": "P12345", "kegg_gene_id": "hsa:1"},
        {"group_id": "up", "input_id": "P12345", "kegg_gene_id": "hsa:1"},
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "up", "input_id": "MISSING"}
    ]


def test_select_groups_resolves_unique_ids_once_then_expands_membership(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))
    mapping_calls = 0
    original_eager_mappings = kegg_module.KEGGDatabase._eager_mappings  # pyright: ignore[reportPrivateUsage]

    def counted_eager_mappings(database: KEGGDatabase) -> pl.DataFrame:
        nonlocal mapping_calls
        mapping_calls += 1
        return original_eager_mappings(database)

    monkeypatch.setattr(
        kegg_module.KEGGDatabase,
        "_eager_mappings",
        counted_eager_mappings,
    )
    selection = db.select_groups(
        {
            " case ": ["sp|P12345|GENE1_HUMAN", "P12345", "MISSING"],
            "control": ["P12345", "MISSING"],
            "empty": [],
        },
        namespace="uniprot",
    )

    groups = selection._df_groups  # pyright: ignore[reportPrivateUsage]
    assert groups is not None
    assert groups["group_id"].to_list() == [
        "case",
        "control",
        "empty",
    ]
    assert selection._df_input_ids["input_id"].to_list() == [  # pyright: ignore[reportPrivateUsage]
        "MISSING",
        "P12345",
    ]
    assert selection.mappings().collect().select(
        "group_id", "input_id", "kegg_pathway_id"
    ).to_dicts() == [
        {
            "group_id": "case",
            "input_id": "P12345",
            "kegg_pathway_id": "hsa00010",
        },
        {
            "group_id": "case",
            "input_id": "P12345",
            "kegg_pathway_id": "hsa01100",
        },
        {
            "group_id": "control",
            "input_id": "P12345",
            "kegg_pathway_id": "hsa00010",
        },
        {
            "group_id": "control",
            "input_id": "P12345",
            "kegg_pathway_id": "hsa01100",
        },
    ]
    selection.mappings().collect()
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "case", "input_id": "MISSING"},
        {"group_id": "control", "input_id": "MISSING"},
    ]
    assert mapping_calls


def test_optional_mapping_files_leave_nullable_columns(tmp_path: Path) -> None:
    files = write_kegg_mapping_fixture(tmp_path)
    db = KEGGDatabase.from_mapping_files(
        uniprot_conversion=files["conv_uniprot"],
        gene_ko=files["gene_ko"],
        gene_pathway=files["gene_pathway"],
        organism_code="hsa",
    )

    row = (
        db.select_ids(["Q9Y243"], namespace="uniprot")
        .mappings()
        .collect()
        .row(
            0,
            named=True,
        )
    )
    assert row["ncbi_gene_id"] is None
    assert row["gene_symbol"] is None
    assert row["gene_description"] is None


def test_write_duckdb_reopens_mapping_without_sidecar(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    path = tmp_path / "kegg.duckdb"
    result = db.write_duckdb(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    reopened = KEGGDatabase.from_duckdb(path)
    assert reopened.mappings().collect().equals(db.mappings().collect())
    with reopened.connect() as connection:
        assert connection.execute("SELECT count(*) FROM mapping").fetchone() == (3,)


def test_mapping_validates_kind_and_snapshot_kind(tmp_path: Path) -> None:
    files = write_kegg_mapping_fixture(tmp_path)
    db = create_mapping_db(files)

    with pytest.raises(ValueError, match="namespace"):
        db.select_ids(["P12345"], namespace="symbol")  # type: ignore[arg-type]

    file_brite = tmp_path / "br08901.json"
    file_brite.write_text('{"name": "ko00001", "children": []}', encoding="utf-8")
    db_brite = KEGGDatabase.from_brite_json(file_brite)
    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="BRITE"):
        db_brite.mappings().collect()


def test_mapping_validates_organism_code(tmp_path: Path) -> None:
    files = write_kegg_mapping_fixture(tmp_path)
    files["gene_ko"].write_text("mmu:1\tko:K00001\n", encoding="utf-8")
    db = create_mapping_db(files)

    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="organism_code"):
        db.mappings().collect()
