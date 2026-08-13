from __future__ import annotations

import os
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.errors import CapabilityError, IntegrityError
from bioextract.reactome import ReactomeDatabase


def write_reactome_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    file_mapping = tmp_path / "UniProt2Reactome.txt"
    file_pathways = tmp_path / "ReactomePathways.txt"
    file_relations = tmp_path / "ReactomePathwaysRelation.txt"

    file_mapping.write_text(
        "\n".join(
            [
                "P04637\tR-HSA-69563\thttps://reactome.org/PathwayBrowser/#/R-HSA-69563\tp53-Dependent G1 DNA Damage Response\tTAS\tHomo sapiens",
                "P04637\tR-HSA-6798695\thttps://reactome.org/PathwayBrowser/#/R-HSA-6798695\tNeutrophil degranulation\tTAS\tHomo sapiens",
                "Q9Y243\tR-HSA-6798695\thttps://reactome.org/PathwayBrowser/#/R-HSA-6798695\tNeutrophil degranulation\tTAS\tHomo sapiens",
                "P31749\tR-MMU-1257604\thttps://reactome.org/PathwayBrowser/#/R-MMU-1257604\tPIP3 activates AKT signaling\tTAS\tMus musculus",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    file_pathways.write_text(
        "\n".join(
            [
                "R-HSA-69563\tp53-Dependent G1 DNA Damage Response\tHomo sapiens",
                "R-HSA-6798695\tNeutrophil degranulation\tHomo sapiens",
                "R-HSA-1640170\tCell Cycle\tHomo sapiens",
                "R-MMU-1257604\tPIP3 activates AKT signaling\tMus musculus",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    file_relations.write_text(
        "\n".join(
            [
                "R-HSA-1640170\tR-HSA-69563",
                "R-HSA-1640170\tR-HSA-6798695",
                "R-MMU-000001\tR-MMU-1257604",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return file_mapping, file_pathways, file_relations


def test_extract_mapping_and_unmapped_single_selection(tmp_path: Path) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)

    selection = (
        ReactomeDatabase.from_files(
            uniprot_mapping=file_mapping,
            pathways=file_pathways,
            relations=file_relations,
        )
        .with_species("Homo sapiens")
        .select_ids([" sp|P04637|P53_HUMAN ", "Q9Y243", "MISSING", ""])
    )

    assert selection.mappings().collect().to_dicts() == [
        {
            "input_id": "P04637",
            "uniprot_id": "P04637",
            "reactome_pathway_id": "R-HSA-6798695",
            "pathway_name": "Neutrophil degranulation",
            "evidence_code": "TAS",
            "species": "Homo sapiens",
            "reactome_url": "https://reactome.org/PathwayBrowser/#/R-HSA-6798695",
        },
        {
            "input_id": "P04637",
            "uniprot_id": "P04637",
            "reactome_pathway_id": "R-HSA-69563",
            "pathway_name": "p53-Dependent G1 DNA Damage Response",
            "evidence_code": "TAS",
            "species": "Homo sapiens",
            "reactome_url": "https://reactome.org/PathwayBrowser/#/R-HSA-69563",
        },
        {
            "input_id": "Q9Y243",
            "uniprot_id": "Q9Y243",
            "reactome_pathway_id": "R-HSA-6798695",
            "pathway_name": "Neutrophil degranulation",
            "evidence_code": "TAS",
            "species": "Homo sapiens",
            "reactome_url": "https://reactome.org/PathwayBrowser/#/R-HSA-6798695",
        },
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [{"input_id": "MISSING"}]


def test_grouped_selection_preserves_groups(tmp_path: Path) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).with_species("Homo sapiens")

    selection = db.select_groups(
        {
            "TumorA": ["P04637", "MISSING"],
            "TumorB": ["P04637", "Q9Y243"],
        }
    )

    df_mapping = selection.mappings().collect()
    assert df_mapping.columns == [
        "group_id",
        "input_id",
        "uniprot_id",
        "reactome_pathway_id",
        "pathway_name",
        "evidence_code",
        "species",
        "reactome_url",
    ]
    assert df_mapping.filter(df_mapping["input_id"] == "P04637").height == 4
    assert selection.mappings().collect().equals(df_mapping)
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "TumorA", "input_id": "MISSING"}
    ]


def test_extract_enrichment_inputs_and_relations_are_species_scoped(
    tmp_path: Path,
) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).with_species("Homo sapiens")

    assert db.pathway_genes().collect().to_dicts() == [
        {"reactome_pathway_id": "R-HSA-6798695", "uniprot_id": "P04637"},
        {"reactome_pathway_id": "R-HSA-6798695", "uniprot_id": "Q9Y243"},
        {"reactome_pathway_id": "R-HSA-69563", "uniprot_id": "P04637"},
    ]
    assert db.pathway_names().collect().to_dicts() == [
        {
            "reactome_pathway_id": "R-HSA-1640170",
            "pathway_name": "Cell Cycle",
            "species": "Homo sapiens",
        },
        {
            "reactome_pathway_id": "R-HSA-6798695",
            "pathway_name": "Neutrophil degranulation",
            "species": "Homo sapiens",
        },
        {
            "reactome_pathway_id": "R-HSA-69563",
            "pathway_name": "p53-Dependent G1 DNA Damage Response",
            "species": "Homo sapiens",
        },
    ]
    assert db.pathway_relations().collect().to_dicts() == [
        {
            "parent_reactome_pathway_id": "R-HSA-1640170",
            "child_reactome_pathway_id": "R-HSA-6798695",
        },
        {
            "parent_reactome_pathway_id": "R-HSA-1640170",
            "child_reactome_pathway_id": "R-HSA-69563",
        },
    ]


def test_build_tidy_writes_duckdb(tmp_path: Path) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    )

    tidy = db.build_tidy()
    assert set(tidy.frames) == {
        "mapping",
        "pathway",
        "relation",
        "term2gene",
        "term2name",
    }
    result = db.write_duckdb(tmp_path / "reactome.duckdb")
    assert result.tables == ("protein_pathway", "pathway", "pathway_relation")
    assert not (tmp_path / "manifest.json").exists()


def test_duckdb_reopen_preserves_domain_selection_and_native_sql(
    tmp_path: Path,
) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    source = ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).with_species("Homo sapiens")
    expected_mapping = (
        source.select_groups(
            {
                "TumorA": ["P04637", "MISSING"],
                "TumorB": ["P04637", "Q9Y243"],
            }
        )
        .mappings()
        .collect()
    )
    expected_unmatched = (
        source.select_groups(
            {
                "TumorA": ["P04637", "MISSING"],
                "TumorB": ["P04637", "Q9Y243"],
            }
        )
        .unmatched_ids()
        .collect()
    )
    expected_term2gene = source.pathway_genes().collect()
    expected_term2name = source.pathway_names().collect()
    expected_relations = source.pathway_relations().collect()
    publication = tmp_path / "reactome.duckdb"
    ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).write_duckdb(publication)

    reopened = ReactomeDatabase.from_duckdb(publication).with_species("Homo sapiens")
    selection = reopened.select_groups(
        {
            "TumorA": ["P04637", "MISSING"],
            "TumorB": ["P04637", "Q9Y243"],
        }
    )
    assert selection.mappings().collect().equals(expected_mapping)
    assert selection.unmatched_ids().collect().equals(expected_unmatched)
    assert reopened.pathway_genes().collect().equals(expected_term2gene)
    assert reopened.pathway_names().collect().equals(expected_term2name)
    assert reopened.pathway_relations().collect().equals(expected_relations)
    assert set(reopened.build_tidy().frames) == set(source.build_tidy().frames)

    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        assert first.execute(
            "SELECT pathway_name FROM pathway WHERE reactome_pathway_id='R-HSA-1640170'"
        ).fetchone() == ("Cell Cycle",)
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_duckdb_reopen_validates_bounded_physical_contract(tmp_path: Path) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    publication = tmp_path / "reactome.duckdb"
    ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).write_duckdb(publication)

    with duckdb.connect(str(publication)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info SET row_count=999999999 "
            "WHERE table_name='protein_pathway'"
        )
    assert (
        ReactomeDatabase.from_duckdb(publication)
        .select_ids(["P04637"])
        .mappings()
        .collect()
        .height
        == 2
    )

    with duckdb.connect(str(publication)) as connection:
        connection.execute("ALTER TABLE pathway DROP COLUMN species")
    with pytest.raises(IntegrityError, match="table schema"):
        ReactomeDatabase.from_duckdb(publication)


def test_duckdb_reopen_rejects_wrong_identity_inventory_and_replacement(
    tmp_path: Path,
) -> None:
    file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)
    wrong = tmp_path / "wrong.duckdb"
    ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).write_duckdb(wrong)
    with duckdb.connect(str(wrong)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='other' "
            "WHERE key='bioextract.resource_name'"
        )
    with pytest.raises(IntegrityError, match="not a bioextract Reactome"):
        ReactomeDatabase.from_duckdb(wrong)

    corrupt = tmp_path / "corrupt.duckdb"
    ReactomeDatabase.from_files(
        uniprot_mapping=file_mapping,
        pathways=file_pathways,
        relations=file_relations,
    ).write_duckdb(corrupt)
    with duckdb.connect(str(corrupt)) as connection:
        connection.execute("CREATE VIEW unexpected AS SELECT 1 AS value")
    with pytest.raises(IntegrityError, match="table/view inventory"):
        ReactomeDatabase.from_duckdb(corrupt)

    current = tmp_path / "current.duckdb"
    replacement = tmp_path / "replacement.duckdb"
    ReactomeDatabase.from_files(uniprot_mapping=file_mapping).write_duckdb(current)
    reopened = ReactomeDatabase.from_duckdb(current)
    cached_selection = reopened.select_ids(["P04637"])
    assert cached_selection.mappings().collect().height == 2
    with pytest.raises(CapabilityError, match="source-file handle"):
        reopened.build_tidy().write_duckdb(tmp_path / "invalid.duckdb")
    with pytest.raises(
        (CapabilityError, pl.exceptions.ComputeError), match="pathway metadata"
    ):
        reopened.pathway_names().collect()
    ReactomeDatabase.from_files(uniprot_mapping=file_mapping).write_duckdb(replacement)
    os.replace(replacement, current)
    with pytest.raises(IntegrityError, match="was replaced"):
        reopened.connect()
    with pytest.raises(
        (IntegrityError, pl.exceptions.ComputeError), match="was replaced"
    ):
        cached_selection.mappings().collect()
    with pytest.raises(IntegrityError, match="was replaced"):
        reopened.build_tidy()


def test_source_handle_rejects_native_connection(tmp_path: Path) -> None:
    file_mapping, _, _ = write_reactome_fixture(tmp_path)
    with pytest.raises(CapabilityError, match="from_duckdb"):
        ReactomeDatabase.from_files(uniprot_mapping=file_mapping).connect()


def test_mapping_only_snapshot_supports_annotation_and_term2gene(
    tmp_path: Path,
) -> None:
    file_mapping, _, _ = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(uniprot_mapping=file_mapping).with_species(
        "Homo sapiens"
    )

    assert db.pathway_genes().collect().height == 3
    assert db.select_ids(["P04637"]).mappings().collect().height == 2

    tidy = db.build_tidy()
    assert set(tidy.frames) == {"mapping", "term2gene"}
    result = db.write_duckdb(tmp_path / "reactome_mapping.duckdb")
    assert result.tables == ("protein_pathway",)

    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="pathways file"):
        db.pathway_names().collect()
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="relations file"
    ):
        db.pathway_relations().collect()


def test_pathway_only_snapshot_supports_term2name(tmp_path: Path) -> None:
    _, file_pathways, _ = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(pathways=file_pathways).with_species(
        "Homo sapiens"
    )

    assert db.pathway_names().collect().height == 3
    assert set(db.build_tidy().frames) == {"pathway", "term2name"}

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="UniProt2Reactome file"
    ):
        db.pathway_genes().collect()


def test_relation_only_snapshot_supports_unscoped_relations(
    tmp_path: Path,
) -> None:
    _, _, file_relations = write_reactome_fixture(tmp_path)
    db = ReactomeDatabase.from_files(relations=file_relations)

    assert db.pathway_relations().collect().height == 3
    assert set(db.build_tidy().frames) == {"relation"}

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError),
        match="species-scoped relation filtering",
    ):
        db.with_species("Homo sapiens").pathway_relations().collect()

    publication = tmp_path / "reactome_relations.duckdb"
    db.write_duckdb(publication)
    with pytest.raises(
        (CapabilityError, pl.exceptions.ComputeError),
        match="species-scoped relation filtering",
    ):
        ReactomeDatabase.from_duckdb(publication).with_species(
            "Homo sapiens"
        ).pathway_relations().collect()


def test_from_files_rejects_missing_files(tmp_path: Path) -> None:
    _file_mapping, file_pathways, file_relations = write_reactome_fixture(tmp_path)

    with pytest.raises(ValueError, match="At least one Reactome input file"):
        ReactomeDatabase.from_files()

    with pytest.raises(FileNotFoundError):
        ReactomeDatabase.from_files(
            uniprot_mapping=tmp_path / "missing.txt",
            pathways=file_pathways,
            relations=file_relations,
        )
