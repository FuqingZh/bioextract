from __future__ import annotations

from pathlib import Path

import duckdb
import pytest

from bioextract.errors import IntegrityError
from bioextract.reactome import ReactomeDatabase
from bioextract.reactome.constant import ENTITY_COLUMN_MAPPING_REASON


def _write_entity_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    complex_pathways = tmp_path / "Complex_2_Pathway_human.txt"
    complex_pathways.write_text(
        "complex\tpathway\ttop_level_pathway\n"
        "R-ALL-1006146\tR-HSA-2\tR-HSA-1\n"
        "R-ALL-1006146\tR-HSA-2\tR-HSA-1\n"
        "R-BAN-1\tR-HSA-3\tR-HSA-1\n",
        encoding="utf-8",
    )
    ewas_pathways = tmp_path / "Ewas2Pathway_human.txt"
    ewas_pathways.write_text(
        "ewas\tpathway\ttop_level_pathway\nR-BAN-5205700\tR-HSA-2\tR-HSA-1\n",
        encoding="utf-8",
    )
    pathways = tmp_path / "ReactomePathways.txt"
    pathways.write_text(
        "R-HSA-1\tTop\tHomo sapiens\n"
        "R-HSA-2\tChild\tHomo sapiens\n"
        "R-HSA-3\tOther child\tHomo sapiens\n",
        encoding="utf-8",
    )
    relations = tmp_path / "ReactomePathwaysRelation.txt"
    relations.write_text("R-HSA-1\tR-HSA-2\nR-HSA-1\tR-HSA-3\n", encoding="utf-8")
    return complex_pathways, ewas_pathways, pathways, relations


def test_entity_relations_preserve_human_scope_and_exact_header_contract(
    tmp_path: Path,
) -> None:
    complex_pathways, ewas_pathways, pathways, relations = _write_entity_fixture(
        tmp_path
    )
    db = ReactomeDatabase.from_files(
        complex_pathways=complex_pathways,
        ewas_pathways=ewas_pathways,
        pathways=pathways,
        relations=relations,
    )
    assert db.complex_pathways().collect().to_dicts() == [
        {
            "reactome_complex_id": "R-ALL-1006146",
            "reactome_pathway_id": "R-HSA-2",
            "top_level_reactome_pathway_id": "R-HSA-1",
        },
        {
            "reactome_complex_id": "R-BAN-1",
            "reactome_pathway_id": "R-HSA-3",
            "top_level_reactome_pathway_id": "R-HSA-1",
        },
    ]
    assert db.ewas_pathways().collect().height == 1
    assert db.with_species("Homo sapiens").complex_pathways().collect().height == 2
    assert db.with_species("Mus musculus").complex_pathways().collect().columns == [
        "reactome_complex_id",
        "reactome_pathway_id",
        "top_level_reactome_pathway_id",
    ]
    assert db.with_species("Mus musculus").complex_pathways().collect().height == 0

    malformed = tmp_path / "malformed.txt"
    malformed.write_text(
        "pathway\tcomplex\ttop_level_pathway\nR-1\tR-2\tR-1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="exact ordered header"):
        ReactomeDatabase.from_files(
            complex_pathways=malformed
        ).complex_pathways().collect()


def test_entity_publication_records_lineage_warnings_and_reopens(
    tmp_path: Path,
) -> None:
    complex_pathways, ewas_pathways, pathways, relations = _write_entity_fixture(
        tmp_path
    )
    publication = tmp_path / "entity.duckdb"
    source = ReactomeDatabase.from_files(
        complex_pathways=complex_pathways,
        ewas_pathways=ewas_pathways,
        pathways=pathways,
        relations=relations,
        release_version="96",
    )
    result = source.write_duckdb(publication)
    assert result.tables == (
        "pathway",
        "pathway_relation",
        "complex_pathway",
        "ewas_pathway",
    )
    with duckdb.connect(str(publication), read_only=True) as connection:
        mappings = connection.execute(
            "SELECT table_name, source_column, output_column, reason "
            "FROM _bioextract.column_mapping ORDER BY table_name, source_column"
        ).fetchall()
        assert mappings == [
            (
                "complex_pathway",
                "complex",
                "reactome_complex_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
            (
                "complex_pathway",
                "pathway",
                "reactome_pathway_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
            (
                "complex_pathway",
                "top_level_pathway",
                "top_level_reactome_pathway_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
            (
                "ewas_pathway",
                "ewas",
                "reactome_ewas_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
            (
                "ewas_pathway",
                "pathway",
                "reactome_pathway_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
            (
                "ewas_pathway",
                "top_level_pathway",
                "top_level_reactome_pathway_id",
                ENTITY_COLUMN_MAPPING_REASON,
            ),
        ]
    reopened = ReactomeDatabase.from_duckdb(publication)
    assert (
        reopened.complex_pathways()
        .collect()
        .equals(source.complex_pathways().collect())
    )
    assert reopened.ewas_pathways().collect().equals(source.ewas_pathways().collect())
    assert reopened.release_version == "96"


def test_entity_top_level_mismatch_is_fatal_when_hierarchy_is_available(
    tmp_path: Path,
) -> None:
    complex_pathways, _, pathways, relations = _write_entity_fixture(tmp_path)
    complex_pathways.write_text(
        "complex\tpathway\ttop_level_pathway\nR-BAN-1\tR-HSA-2\tR-HSA-3\n",
        encoding="utf-8",
    )
    db = ReactomeDatabase.from_files(
        complex_pathways=complex_pathways,
        pathways=pathways,
        relations=relations,
    )
    with pytest.raises(IntegrityError, match="not an ancestor"):
        db.write_duckdb(tmp_path / "invalid.duckdb")
