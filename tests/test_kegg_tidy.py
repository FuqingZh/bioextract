import json
from pathlib import Path

import polars as pl

from bioextract.kegg import KeggDb
from bioextract.kegg.brite import run_tidy_kegg_brite
from bioextract.kegg.brite.parse import parse_entry_and_ko, parse_pathway_level3


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


def test_parse_pathway_level3_supports_path_br_and_plain_forms() -> None:
    level3_path = parse_pathway_level3(
        "00010 Glycolysis / Gluconeogenesis [PATH:tcar00010]"
    )
    assert level3_path.id == "00010"
    assert level3_path.kegg_id == "tcar00010"
    assert level3_path.name == "Glycolysis / Gluconeogenesis"

    level3_plain = parse_pathway_level3("00566 Sulfoquinovose metabolism")
    assert level3_plain.id == "00566"
    assert level3_plain.kegg_id is None
    assert level3_plain.name == "Sulfoquinovose metabolism"

    level3_br = parse_pathway_level3("01001 Protein kinases [BR:tcar01001]")
    assert level3_br.id == "01001"
    assert level3_br.kegg_id == "tcar01001"
    assert level3_br.name == "Protein kinases"


def test_parse_entry_and_ko_supports_entry_only_leaf() -> None:
    leaf = parse_entry_and_ko("U0034_01605 oxidoreductase")
    assert leaf.entry.id == "U0034_01605"
    assert leaf.entry.name == "oxidoreductase"
    assert leaf.ko is None


def test_kegg_db_build_tidy_exposes_frames_and_write_contract(tmp_path: Path) -> None:
    file_in = tmp_path / "tcar00001.json"
    dir_out = tmp_path / "tidy"
    write_minimal_brite_json(file_in)

    tidy = KeggDb.from_brite_json(file_in).build_tidy()

    assert set(tidy.frames) == {"pathway"}
    assert tidy.frames["pathway"].height == 2

    report = tidy.write(dir_out)

    assert report.manifest is None
    assert len(report.assets) == 1
    assert (dir_out / "pathway.parquet").exists()
    assert not (dir_out / "manifest.json").exists()

    report_manifest = tidy.write(dir_out / "with_manifest", should_write_manifest=True)
    assert report_manifest.manifest is not None
    assert report_manifest.manifest["schema_version"] == "kegg-brite-tidy-v0.1"
    data_manifest = json.loads(
        (dir_out / "with_manifest" / "manifest.json").read_text("utf-8")
    )
    assert data_manifest["sources"][0]["path"] == file_in.as_posix()
    assert data_manifest["sources"][0]["media_type"] == "application/json"

    df_pathway = pl.read_parquet(dir_out / "pathway.parquet")
    assert df_pathway.to_dicts()[0] == {
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


def test_legacy_kegg_tidy_runner_still_writes_contract(tmp_path: Path) -> None:
    file_in = tmp_path / "tcar00001.json"
    dir_out = tmp_path / "legacy"
    write_minimal_brite_json(file_in, has_leaf=False)

    run_tidy_kegg_brite(file_in=file_in, dir_out=dir_out)

    df_pathway = pl.read_parquet(dir_out / "pathway.parquet")
    assert df_pathway.height == 1
    assert df_pathway.to_dicts()[0]["entry_id"] is None
