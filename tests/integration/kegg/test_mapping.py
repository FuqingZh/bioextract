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

    df_mapping = db.extract_mapping()

    assert df_mapping.columns == [
        "OrganismCode",
        "KeggGeneId",
        "UniProtId",
        "NcbiGeneId",
        "KoId",
        "KeggPathwayId",
        "PathwayMapId",
        "GeneSymbol",
        "GeneDescription",
    ]
    assert df_mapping.to_dicts() == [
        {
            "OrganismCode": "hsa",
            "KeggGeneId": "hsa:1",
            "UniProtId": "P12345",
            "NcbiGeneId": "101",
            "KoId": "K00001",
            "KeggPathwayId": "hsa00010",
            "PathwayMapId": "map00010",
            "GeneSymbol": "GENE1",
            "GeneDescription": "alpha description",
        },
        {
            "OrganismCode": "hsa",
            "KeggGeneId": "hsa:1",
            "UniProtId": "P12345",
            "NcbiGeneId": "101",
            "KoId": "K00001",
            "KeggPathwayId": "hsa01100",
            "PathwayMapId": "map01100",
            "GeneSymbol": "GENE1",
            "GeneDescription": "alpha description",
        },
        {
            "OrganismCode": "hsa",
            "KeggGeneId": "hsa:2",
            "UniProtId": "Q9Y243",
            "NcbiGeneId": "102",
            "KoId": "K00002",
            "KeggPathwayId": "hsa04110",
            "PathwayMapId": "map04110",
            "GeneSymbol": "GENE2",
            "GeneDescription": None,
        },
    ]


def test_select_ids_supports_input_id_kinds_and_unmapped(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    df_uniprot = db.select_ids(
        ["sp|P12345|GENE1_HUMAN", "MISSING"],
        namespace="uniprot",
    ).extract_mapping()
    assert df_uniprot.select("InputId", "InputNamespace", "KeggGeneId").to_dicts() == [
        {"InputId": "P12345", "InputNamespace": "uniprot", "KeggGeneId": "hsa:1"},
        {"InputId": "P12345", "InputNamespace": "uniprot", "KeggGeneId": "hsa:1"},
    ]
    assert db.select_ids(
        ["P12345", "MISSING"],
        namespace="uniprot",
    ).extract_unmatched_ids().to_dicts() == [{"InputId": "MISSING"}]

    df_ncbi = db.select_ids(["102"], namespace="ncbi_gene").extract_mapping()
    assert df_ncbi.select("InputId", "KeggGeneId").to_dicts() == [
        {"InputId": "102", "KeggGeneId": "hsa:2"}
    ]

    df_kegg = db.select_ids(["hsa:1"], namespace="kegg_gene").extract_mapping()
    assert df_kegg.select("InputId", "UniProtId", "KeggPathwayId").to_dicts() == [
        {"InputId": "hsa:1", "UniProtId": "P12345", "KeggPathwayId": "hsa00010"},
        {"InputId": "hsa:1", "UniProtId": "P12345", "KeggPathwayId": "hsa01100"},
    ]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    selection = db.select_groups(
        {"up": ["P12345", "MISSING"], "down": ["Q9Y243"]},
        namespace="uniprot",
    )

    df_mapping = selection.extract_mapping()
    assert df_mapping.columns[:3] == ["GroupId", "InputId", "InputNamespace"]
    assert df_mapping.select("GroupId", "InputId", "KeggGeneId").to_dicts() == [
        {"GroupId": "down", "InputId": "Q9Y243", "KeggGeneId": "hsa:2"},
        {"GroupId": "up", "InputId": "P12345", "KeggGeneId": "hsa:1"},
        {"GroupId": "up", "InputId": "P12345", "KeggGeneId": "hsa:1"},
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"GroupId": "up", "InputId": "MISSING"}
    ]


def test_select_groups_resolves_unique_ids_once_then_expands_membership(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))
    mapping_calls = 0
    original_extract_mapping = kegg_module.KEGGDatabase.extract_mapping

    def counted_extract_mapping(database: KEGGDatabase) -> pl.DataFrame:
        nonlocal mapping_calls
        mapping_calls += 1
        return original_extract_mapping(database)

    monkeypatch.setattr(
        kegg_module.KEGGDatabase,
        "extract_mapping",
        counted_extract_mapping,
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
    assert groups["GroupId"].to_list() == [
        "case",
        "control",
        "empty",
    ]
    assert selection._df_input_ids["InputId"].to_list() == [  # pyright: ignore[reportPrivateUsage]
        "MISSING",
        "P12345",
    ]
    assert selection.extract_mapping().select(
        "GroupId", "InputId", "KeggPathwayId"
    ).to_dicts() == [
        {
            "GroupId": "case",
            "InputId": "P12345",
            "KeggPathwayId": "hsa00010",
        },
        {
            "GroupId": "case",
            "InputId": "P12345",
            "KeggPathwayId": "hsa01100",
        },
        {
            "GroupId": "control",
            "InputId": "P12345",
            "KeggPathwayId": "hsa00010",
        },
        {
            "GroupId": "control",
            "InputId": "P12345",
            "KeggPathwayId": "hsa01100",
        },
    ]
    selection.extract_mapping()
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"GroupId": "case", "InputId": "MISSING"},
        {"GroupId": "control", "InputId": "MISSING"},
    ]
    assert mapping_calls == 1


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
        .extract_mapping()
        .row(
            0,
            named=True,
        )
    )
    assert row["NcbiGeneId"] is None
    assert row["GeneSymbol"] is None
    assert row["GeneDescription"] is None


def test_write_duckdb_reopens_mapping_without_sidecar(tmp_path: Path) -> None:
    db = create_mapping_db(write_kegg_mapping_fixture(tmp_path))

    path = tmp_path / "kegg.duckdb"
    result = db.write_duckdb(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    reopened = KEGGDatabase.from_duckdb(path)
    assert reopened.extract_mapping().equals(db.extract_mapping())
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
    with pytest.raises(ValueError, match="BRITE JSON snapshot"):
        db_brite.extract_mapping()


def test_mapping_validates_organism_code(tmp_path: Path) -> None:
    files = write_kegg_mapping_fixture(tmp_path)
    files["gene_ko"].write_text("mmu:1\tko:K00001\n", encoding="utf-8")
    db = create_mapping_db(files)

    with pytest.raises(ValueError, match="organism_code"):
        db.extract_mapping()
