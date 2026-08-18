from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.reactome import ReactomeDatabase
from bioextract.reactome.constant import (
    MAPPING_ROLE_SPECS,
    PATHWAY_ROLE,
    RELATION_ROLE,
)


def _write_mapping_matrix_fixture(tmp_path: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    source_values = {
        "uniprot": "P04637",
        "ncbi": "NP_001",
        "chebi": "100241",
        "gtop": "1",
    }
    for spec in MAPPING_ROLE_SPECS:
        source_value = source_values[spec.namespace]
        event_value = "R-HSA-1" if spec.target == "pathway" else "R-HSA-2"
        name_value = (
            "Fixture pathway" if spec.target == "pathway" else 'Fixture "reaction"'
        )
        row = [
            source_value,
            event_value,
            f"https://reactome.org/PathwayBrowser/#/{event_value}",
            name_value,
            "TAS",
            "Homo sapiens",
        ]
        path = tmp_path / spec.filename
        path.write_text("\t".join(row) + "\n", encoding="utf-8")
        paths[spec.argument_name] = path

    pathways = tmp_path / "ReactomePathways.txt"
    pathways.write_text("R-HSA-1\tFixture pathway\tHomo sapiens\n", encoding="utf-8")
    paths["pathways"] = pathways
    relations = tmp_path / "ReactomePathwaysRelation.txt"
    relations.write_text("R-HSA-1\tR-HSA-1\n", encoding="utf-8")
    paths["relations"] = relations
    return paths


def test_mapping_role_registry_is_complete_and_private() -> None:
    assert len(MAPPING_ROLE_SPECS) == 12
    assert len({spec.role for spec in MAPPING_ROLE_SPECS}) == 12
    assert len({spec.argument_name for spec in MAPPING_ROLE_SPECS}) == 12
    assert len({spec.filename for spec in MAPPING_ROLE_SPECS}) == 12
    assert len({spec.source_column for spec in MAPPING_ROLE_SPECS}) == 4
    assert {spec.namespace for spec in MAPPING_ROLE_SPECS} == {
        "uniprot",
        "ncbi",
        "chebi",
        "gtop",
    }


def test_all_matrix_roles_have_distinct_whole_resource_relations(
    tmp_path: Path,
) -> None:
    paths = _write_mapping_matrix_fixture(tmp_path)
    db = ReactomeDatabase.from_files(**paths)  # pyright: ignore[reportArgumentType]

    assert set(db.build_tidy().frames) == {spec.role for spec in MAPPING_ROLE_SPECS} | {
        PATHWAY_ROLE,
        RELATION_ROLE,
    }
    for spec in MAPPING_ROLE_SPECS:
        if spec.target == "pathway":
            assert spec.pathway_level is not None
            frame = db.pathway_mappings(
                namespace=spec.namespace,
                pathway_level=spec.pathway_level,
            ).collect()
        else:
            frame = db.reaction_mappings(namespace=spec.namespace).collect()
        assert frame.columns == list(spec.public_columns)
        assert frame.height == 1

    publication = tmp_path / "matrix.duckdb"
    result = db.write_duckdb(publication)
    assert result.tables == tuple(
        [spec.role for spec in MAPPING_ROLE_SPECS] + [PATHWAY_ROLE, RELATION_ROLE]
    )
    reopened = ReactomeDatabase.from_duckdb(publication)
    assert reopened.reaction_mappings(namespace="chebi").collect().height == 1
    assert (
        reopened.pathway_mappings(namespace="ncbi", pathway_level="all_levels")
        .collect()
        .height
        == 1
    )
    with duckdb.connect(str(publication), read_only=True) as connection:
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key='bioextract.resource_schema_version'"
        ).fetchone() == ("reactome-mapping-v0.5",)


def test_namespace_specific_selection_lineage_and_reaction_shape(
    tmp_path: Path,
) -> None:
    paths = _write_mapping_matrix_fixture(tmp_path)
    db = ReactomeDatabase.from_files(  # pyright: ignore[reportArgumentType]
        **paths  # pyright: ignore[reportArgumentType]
    ).with_species("Homo sapiens")

    ncbi = db.select_ids([" NP_001 "], namespace="ncbi")
    assert ncbi.mappings().collect().to_dicts() == [
        {
            "input_id": "NP_001",
            "ncbi_id": "NP_001",
            "reactome_pathway_id": "R-HSA-1",
            "pathway_name": "Fixture pathway",
            "evidence_code": "TAS",
            "species": "Homo sapiens",
            "reactome_url": "https://reactome.org/PathwayBrowser/#/R-HSA-1",
        }
    ]

    chebi = db.select_ids(["100241", "CHEBI:100241"], namespace="chebi")
    assert chebi.mappings().collect().get_column("input_id").to_list() == [
        "CHEBI:100241"
    ]

    reaction = db.select_groups(
        {"case": ["CHEBI:100241"], "control": ["CHEBI:999999"]},
        namespace="chebi",
        target="reaction",
    )
    assert reaction.mappings().collect().columns == [
        "group_id",
        "input_id",
        "chebi_id",
        "reactome_reaction_id",
        "reaction_name",
        "evidence_code",
        "species",
        "reactome_url",
    ]
    assert reaction.unmatched_ids().collect().to_dicts() == [
        {"group_id": "control", "input_id": "CHEBI:999999"}
    ]


@pytest.mark.parametrize(
    "selection_kwargs",
    [
        pytest.param(
            {"target": "pathway", "pathway_level": "lowest_level"},
            id="pathway-lowest-level",
        ),
        pytest.param(
            {"target": "pathway", "pathway_level": "all_levels"},
            id="pathway-all-levels",
        ),
        pytest.param({"target": "reaction"}, id="reaction"),
    ],
)
@pytest.mark.parametrize(
    ("namespace", "matched_id", "missing_id"),
    [
        ("uniprot", "P04637", "P999999"),
        ("ncbi", "NP_001", "NP_missing"),
        ("chebi", "CHEBI:100241", "CHEBI:999999"),
        ("gtop", "1", "999999"),
    ],
)
def test_publication_single_selection_unmatched_ids_keep_public_schema(
    tmp_path: Path,
    selection_kwargs: dict[str, str],
    namespace: str,
    matched_id: str,
    missing_id: str,
) -> None:
    paths = _write_mapping_matrix_fixture(tmp_path)
    source = ReactomeDatabase.from_files(  # pyright: ignore[reportArgumentType]
        **paths  # pyright: ignore[reportArgumentType]
    ).with_species("Homo sapiens")
    publication = tmp_path / "matrix-selection.duckdb"
    source.write_duckdb(publication)
    reopened = ReactomeDatabase.from_duckdb(publication).with_species("Homo sapiens")

    for database in (source, reopened):
        selection = database.select_ids(
            [matched_id, missing_id],
            namespace=namespace,
            **selection_kwargs,
        )
        assert selection.mappings().collect().height == 1
        unmatched = selection.unmatched_ids().collect()
        assert unmatched.schema == {"input_id": pl.String}
        assert unmatched.to_dicts() == [{"input_id": missing_id}]


def test_reaction_pathway_level_and_namespace_grammars_fail_closed(
    tmp_path: Path,
) -> None:
    paths = _write_mapping_matrix_fixture(tmp_path)
    db = ReactomeDatabase.from_files(**paths)  # pyright: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="pathway_level"):
        db.select_ids(
            ["1"], namespace="gtop", target="reaction", pathway_level="lowest_level"
        )
    with pytest.raises(ValueError, match="Invalid GtoP"):
        db.select_ids(["GtoP:1"], namespace="gtop")
    with pytest.raises(ValueError, match="Invalid ChEBI"):
        db.select_ids(["CHEBI:abc"], namespace="chebi")
