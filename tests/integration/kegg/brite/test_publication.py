import json
from pathlib import Path

import polars as pl

from bioextract.kegg import KEGGDatabase


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
