import json
from pathlib import Path

import polars as pl

from bioextract.go import GoDb
from bioextract.go.ontology import run_tidy_go_ontology


def write_minimal_obo(file_in: Path) -> None:
    file_in.write_text(
        """format-version: 1.2

[Term]
id: GO:0000001
name: root process
namespace: biological_process
def: "root definition" [GOC:ai]
comment: root comment
xref: Wikipedia:Root
synonym: "root proc" EXACT []

[Term]
id: GO:0000002
name: child process
namespace: biological_process
def: "child definition" [GOC:ai]
is_a: GO:0000001 ! root process
relationship: part_of GO:0000001 ! root process
alt_id: GO:1234567
xref: Reactome:R-HSA-12345
synonym: "child proc" RELATED [GOC:ai]

[Term]
id: GO:0000003
name: leaf function
namespace: molecular_function
def: "leaf definition" [GOC:ai]
is_a: GO:0000002 ! child process
relationship: regulates GO:0000001 ! root process

[Term]
id: GO:0000004
name: obsolete component
namespace: cellular_component
def: "obsolete definition" [GOC:ai]
is_obsolete: true

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component
def: "The part of a cell or its extracellular environment." [GOC:go_curators]

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
def: "All of the contents of a cell excluding the plasma membrane and nucleus." [GOC:go_curators]
is_a: GO:0005575 ! cellular_component
""",
        encoding="utf-8",
    )


def test_go_db_build_tidy_exposes_frames_and_write_contract(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    dir_out = tmp_path / "tidy"
    write_minimal_obo(file_in)

    tidy = GoDb.from_obo(file_in).build_tidy()

    assert set(tidy.frames) == {
        "term",
        "edge",
        "synonym",
        "xref",
        "alt_id",
        "ancestor_all",
        "depth",
    }
    assert tidy.frames["term"].height == 6
    assert tidy.frames["edge"].height == 5

    report = tidy.write(dir_out)

    assert report.dir_out == dir_out
    assert report.manifest is None
    assert len(report.assets) == 7
    assert (dir_out / "term.parquet").exists()
    assert (dir_out / "ancestor_all.parquet").exists()
    assert not (dir_out / "manifest.json").exists()

    report_manifest = tidy.write(dir_out / "with_manifest", should_write_manifest=True)
    assert report_manifest.manifest is not None
    assert report_manifest.manifest["schema_version"] == "go-obo-tidy-v0.1"
    data_manifest = json.loads(
        (dir_out / "with_manifest" / "manifest.json").read_text("utf-8")
    )
    assert data_manifest["sources"][0]["path"] == file_in.as_posix()
    assert data_manifest["sources"][0]["media_type"] == "text/obo"

    df_term = pl.read_parquet(dir_out / "term.parquet")
    row_child = (
        df_term.filter(pl.col("go_id") == "GO:0000002")
        .select("term_name", "definition", "is_obsolete")
        .to_dicts()[0]
    )
    assert row_child == {
        "term_name": "child process",
        "definition": "child definition",
        "is_obsolete": False,
    }


def test_legacy_go_tidy_runner_still_writes_contract(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    dir_out = tmp_path / "legacy"
    write_minimal_obo(file_in)

    run_tidy_go_ontology(file_in=file_in, dir_out=dir_out)

    assert not (dir_out / "manifest.json").exists()
    assert pl.read_parquet(dir_out / "term.parquet").height == 6


def test_go_db_extracts_subcell_from_cellular_component(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    df_subcell = GoDb.from_obo(file_in).extract_subcell()

    assert df_subcell.to_dicts() == [
        {
            "go_id": "GO:0005575",
            "subcell_name": "cellular_component",
            "definition": "The part of a cell or its extracellular environment.",
            "min_depth_from_root": 0,
            "max_depth_from_root": 0,
        },
        {
            "go_id": "GO:0005737",
            "subcell_name": "cytoplasm",
            "definition": (
                "All of the contents of a cell excluding the plasma membrane and "
                "nucleus."
            ),
            "min_depth_from_root": 1,
            "max_depth_from_root": 1,
        },
    ]

    df_subcell_with_obsolete = GoDb.from_obo(file_in).extract_subcell(
        include_obsolete=True
    )
    assert df_subcell_with_obsolete["go_id"].to_list() == [
        "GO:0000004",
        "GO:0005575",
        "GO:0005737",
    ]


def test_go_db_writes_subcell_parquet(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    file_out = tmp_path / "subcell.parquet"
    write_minimal_obo(file_in)

    path_written = GoDb.from_obo(file_in).write_subcell(file_out)

    assert path_written == file_out
    assert pl.read_parquet(file_out).height == 2
