from __future__ import annotations

from pathlib import Path

import pytest

from bioextract.wikipathways import WikiPathwaysDatabase


def write_wikipathways_fixture(tmp_path: Path) -> Path:
    file_gmt = tmp_path / "wikipathways-20260510-gmt-Homo_sapiens.gmt"
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


def test_extract_pathway_term_frames_and_species_filter(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt, species="Homo sapiens")

    assert db.extract_pathway().to_dicts() == [
        {
            "WikiPathwaysId": "WP100",
            "PathwayName": "Glutathione metabolism",
            "Species": "Homo sapiens",
            "Collection": "WikiPathways_20260510",
            "Version": "20260510",
            "Url": "https://www.wikipathways.org/instance/WP100",
            "GeneCount": 2,
        },
        {
            "WikiPathwaysId": "WP106",
            "PathwayName": "Alanine and aspartate metabolism",
            "Species": "Homo sapiens",
            "Collection": "WikiPathways_20260510",
            "Version": "20260510",
            "Url": "https://www.wikipathways.org/instance/WP106",
            "GeneCount": 2,
        },
    ]
    assert db.extract_term2gene().to_dicts() == [
        {"WikiPathwaysId": "WP100", "GeneId": "2678"},
        {"WikiPathwaysId": "WP100", "GeneId": "2687"},
        {"WikiPathwaysId": "WP106", "GeneId": "2806"},
        {"WikiPathwaysId": "WP106", "GeneId": "435"},
    ]
    assert db.extract_term2name().columns == [
        "WikiPathwaysId",
        "PathwayName",
        "Species",
        "Collection",
        "Version",
        "Url",
    ]


def test_single_and_grouped_selection(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt, species="Homo sapiens")

    selection = db.select_ids(["2687", " 435 ", "MISSING", ""])
    assert selection.extract_mapping().to_dicts() == [
        {
            "InputId": "2687",
            "GeneId": "2687",
            "WikiPathwaysId": "WP100",
            "PathwayName": "Glutathione metabolism",
            "Species": "Homo sapiens",
            "Url": "https://www.wikipathways.org/instance/WP100",
        },
        {
            "InputId": "435",
            "GeneId": "435",
            "WikiPathwaysId": "WP106",
            "PathwayName": "Alanine and aspartate metabolism",
            "Species": "Homo sapiens",
            "Url": "https://www.wikipathways.org/instance/WP106",
        },
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [{"InputId": "MISSING"}]

    grouped = db.select_groups({"A": ["2687"], "B": ["2687", "MISSING"]})
    assert grouped.extract_mapping().columns[0] == "GroupId"
    assert grouped.extract_mapping().height == 2
    assert grouped.extract_unmatched_ids().to_dicts() == [
        {"GroupId": "B", "InputId": "MISSING"}
    ]


def test_build_tidy_writes_duckdb_without_sidecar(tmp_path: Path) -> None:
    file_gmt = write_wikipathways_fixture(tmp_path)
    db = WikiPathwaysDatabase.from_gmt(file_gmt, species="Homo sapiens")

    tidy = db.build_tidy()
    assert set(tidy.frames) == {"pathway", "term2gene", "term2name"}
    result = db.write_duckdb(tmp_path / "wikipathways.duckdb")
    assert result.tables == ("pathway", "pathway_gene")
    assert not (tmp_path / "manifest.json").exists()


def test_from_gmt_rejects_missing_and_malformed_files(
    tmp_path: Path,
) -> None:
    write_wikipathways_fixture(tmp_path)

    with pytest.raises(FileNotFoundError):
        WikiPathwaysDatabase.from_gmt(tmp_path / "missing.gmt")

    file_bad = tmp_path / "bad.gmt"
    file_bad.write_text("bad\thttps://example.org\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="four '%' separated fields"):
        WikiPathwaysDatabase.from_gmt(file_bad).extract_pathway()
