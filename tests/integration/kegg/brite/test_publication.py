import json
from importlib import import_module
from pathlib import Path
from typing import Protocol, cast

import polars as pl
import pytest

from bioextract.errors import IntegrityError
from bioextract.kegg import KEGGDatabase


class _TidyValidator(Protocol):
    def __call__(self, path: Path, *, profile: str) -> object: ...


def write_minimal_brite_json(file_in: Path, *, has_leaf: bool = True) -> None:
    leaf_nodes = (
        [
            {
                "name": (
                    "U0034_04525 "
                    "bifunctional transcriptional regulator/glucokinase\t"
                    "K00845 glk; glucokinase [EC:2.7.1.2]"
                )
            },
            {
                "name": (
                    "U0034_00675 "
                    "pgi; glucose-6-phosphate isomerase\t"
                    "K01810 GPI; glucose-6-phosphate isomerase [EC:5.3.1.9]"
                )
            },
        ]
        if has_leaf
        else []
    )
    level3_name = (
        "00010 Glycolysis / Gluconeogenesis [PATH:tcar00010]"
        if has_leaf
        else "00566 Sulfoquinovose metabolism"
    )
    file_in.write_text(
        json.dumps(
            {
                "name": "tcar00001",
                "children": [
                    {
                        "name": "09100 Metabolism",
                        "children": [
                            {
                                "name": "09101 Carbohydrate metabolism",
                                "children": [
                                    {
                                        "name": level3_name,
                                        "children": leaf_nodes,
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


def test_kegg_db_build_tidy_exposes_frames_and_writes_duckdb(tmp_path: Path) -> None:
    file_in = tmp_path / "tcar00001.json"
    path = tmp_path / "kegg.duckdb"
    write_minimal_brite_json(file_in)

    db = KEGGDatabase.from_brite_json(file_in)
    tidy = db.build_tidy()

    assert set(tidy.frames) == {"pathway"}
    assert tidy.frames["pathway"].select(pl.len()).collect().item() == 2

    result = db.write_duckdb(path)
    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    reopened = KEGGDatabase.from_duckdb(path)
    with reopened.connect() as connection:
        columns = connection.sql("FROM pathway").columns
        rows = connection.execute(
            "SELECT * FROM pathway WHERE entry_id='U0034_04525'"
        ).fetchall()
    assert columns == list(tidy.frames["pathway"].collect_schema())
    assert dict(zip(columns, rows[0], strict=True)) == {
        "pathway_level1_id": "09100",
        "pathway_level1_name": "Metabolism",
        "pathway_level2_id": "09101",
        "pathway_level2_name": "Carbohydrate metabolism",
        "pathway_level3_id": "00010",
        "pathway_level3_kegg_id": "tcar00010",
        "pathway_level3_name": "Glycolysis / Gluconeogenesis",
        "entry_id": "U0034_04525",
        "entry_name": "bifunctional transcriptional regulator/glucokinase",
        "ko_id": "K00845",
        "ko_name": "glk; glucokinase [EC:2.7.1.2]",
    }


def test_reopened_brite_handle_rejects_replaced_publication(tmp_path: Path) -> None:
    first_source = tmp_path / "first.json"
    second_source = tmp_path / "second.json"
    target = tmp_path / "kegg.duckdb"
    replacement = tmp_path / "replacement.duckdb"
    write_minimal_brite_json(first_source)
    write_minimal_brite_json(second_source)
    KEGGDatabase.from_brite_json(first_source).write_duckdb(target)
    stale = KEGGDatabase.from_duckdb(target)
    KEGGDatabase.from_brite_json(second_source).write_duckdb(replacement)

    replacement.replace(target)

    with pytest.raises(IntegrityError, match="BRITE publication was replaced"):
        stale.connect()
    fresh = KEGGDatabase.from_duckdb(target)
    with fresh.connect() as connection:
        assert connection.execute("SELECT count(*) FROM pathway").fetchone() == (2,)


def test_from_duckdb_rejects_replacement_during_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_source = tmp_path / "first.json"
    second_source = tmp_path / "second.json"
    target = tmp_path / "kegg.duckdb"
    replacement = tmp_path / "replacement.duckdb"
    write_minimal_brite_json(first_source)
    write_minimal_brite_json(second_source)
    KEGGDatabase.from_brite_json(first_source).write_duckdb(target)
    KEGGDatabase.from_brite_json(second_source).write_duckdb(replacement)
    implementation = import_module("bioextract.kegg.kegg")
    original = cast(
        "_TidyValidator",
        implementation.__dict__["_validate_tidy_publication"],
    )

    def replace_after_validation(path: Path, *, profile: str) -> object:
        validated = original(path, profile=profile)
        replacement.replace(path)
        return validated

    monkeypatch.setattr(
        implementation,
        "_validate_tidy_publication",
        replace_after_validation,
    )

    with pytest.raises(IntegrityError, match="changed during validation"):
        KEGGDatabase.from_duckdb(target)


def test_from_duckdb_classifies_interrupted_replacement_as_integrity_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_source = tmp_path / "first.json"
    second_source = tmp_path / "second.json"
    target = tmp_path / "kegg.duckdb"
    replacement = tmp_path / "replacement.duckdb"
    write_minimal_brite_json(first_source)
    write_minimal_brite_json(second_source)
    KEGGDatabase.from_brite_json(first_source).write_duckdb(target)
    KEGGDatabase.from_brite_json(second_source).write_duckdb(replacement)
    implementation = import_module("bioextract.kegg.kegg")

    def interrupt_validation(path: Path, *, profile: str) -> object:
        del path, profile
        replacement.replace(target)
        raise ValueError("simulated validator interruption")

    monkeypatch.setattr(
        implementation,
        "_validate_tidy_publication",
        interrupt_validation,
    )

    with pytest.raises(IntegrityError, match="changed during validation"):
        KEGGDatabase.from_duckdb(target)
