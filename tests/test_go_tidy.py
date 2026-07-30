import gzip
import tarfile
import zipfile
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.go import GODatabase, GoSubsetId


def write_minimal_obo(file_in: Path) -> None:
    file_in.write_text(
        """format-version: 1.2
subsetdef: goslim_generic "Generic GO slim"
subsetdef: goslim_plant "Plant GO slim"

[Term]
id: GO:0000001
name: root process
namespace: biological_process
def: "root definition" [GOC:ai]
comment: root comment
subset: goslim_generic
xref: Wikipedia:Root
synonym: "root proc" EXACT []

[Term]
id: GO:0000002
name: child process
namespace: biological_process
def: "child definition" [GOC:ai]
subset: goslim_generic
subset: goslim_plant
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
subset: goslim_plant
is_a: GO:0000002 ! child process
relationship: regulates GO:0000001 ! root process

[Term]
id: GO:0000004
name: obsolete component
namespace: cellular_component
def: "obsolete definition" [GOC:ai]
subset: goslim_generic
is_obsolete: true

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component
def: "The part of a cell or its extracellular environment." [GOC:go_curators]
subset: goslim_generic

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
def: "All of the contents of a cell excluding the plasma membrane and nucleus." [GOC:go_curators]
subset: goslim_generic
is_a: GO:0005575 ! cellular_component
""",
        encoding="utf-8",
    )


@pytest.mark.parametrize("container", ["plain", "gzip", "zip", "tar", "directory"])
def test_go_obo_container_is_detected_internally(
    tmp_path: Path,
    container: str,
) -> None:
    file_plain = tmp_path / "go.obo"
    write_minimal_obo(file_plain)
    source = tmp_path / f"go-{container}.snapshot"
    if container == "plain":
        source = file_plain
    elif container == "gzip":
        with (
            file_plain.open("rb") as handle_in,
            gzip.open(source, "wb") as handle_out,
        ):
            handle_out.write(handle_in.read())
    elif container == "zip":
        with zipfile.ZipFile(source, "w") as archive:
            archive.write(file_plain, arcname="ontology/go.obo")
    elif container == "tar":
        with tarfile.open(source, "w") as archive:
            archive.add(file_plain, arcname="ontology/go.obo")
    else:
        source = tmp_path / "release"
        nested = source / "ontology"
        nested.mkdir(parents=True)
        file_plain.replace(nested / "go.obo")

    tidy = GODatabase.from_obo(source).build_tidy()
    assert tidy.frames["term"].select(pl.len()).collect().item() == 6


def test_go_db_build_tidy_exposes_frames_and_writes_duckdb(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    path = tmp_path / "go.duckdb"
    write_minimal_obo(file_in)

    tidy = GODatabase.from_obo(file_in).build_tidy()

    assert set(tidy.frames) == {
        "term",
        "edge",
        "synonym",
        "xref",
        "alt_id",
        "subset_membership",
        "subset_definition",
        "ancestor_all",
        "depth",
    }
    assert tidy.frames["term"].select(pl.len()).collect().item() == 6
    assert tidy.frames["edge"].select(pl.len()).collect().item() == 5
    assert tidy.frames["subset_membership"].select(pl.len()).collect().item() == 7
    assert tidy.frames["subset_definition"].select(pl.len()).collect().item() == 2

    result = GODatabase.from_obo(file_in).write_duckdb(path)
    assert "term_relation" in result.tables
    assert not (tmp_path / "manifest.json").exists()
    with duckdb.connect(str(path), read_only=True) as connection:
        df_term = pl.read_database("SELECT * FROM term", connection)
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


def test_go_db_lists_subsets(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    df_subsets = GODatabase.from_obo(file_in).list_subsets()

    assert df_subsets.to_dicts() == [
        {
            "subset_id": "goslim_generic",
            "subset_name": "Generic GO slim",
            "num_terms": 5,
        },
        {
            "subset_id": "goslim_plant",
            "subset_name": "Plant GO slim",
            "num_terms": 2,
        },
    ]


def test_go_db_selects_terms_by_namespace_subset_and_alt_id(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)
    db = GODatabase.from_obo(file_in)

    df_cellular_generic = db.select_terms(
        namespace="cellular_component",
        subset_id=GoSubsetId.GOSLIM_GENERIC,
    )
    assert df_cellular_generic.select("go_id", "subset_id").to_dicts() == [
        {"go_id": "GO:0005575", "subset_id": "goslim_generic"},
        {"go_id": "GO:0005737", "subset_id": "goslim_generic"},
    ]

    df_selected = db.select_terms(term_ids=["GO:1234567", "GO:0000003"])
    assert df_selected.select("input_go_id", "go_id", "term_name").to_dicts() == [
        {
            "input_go_id": "GO:1234567",
            "go_id": "GO:0000002",
            "term_name": "child process",
        },
        {
            "input_go_id": "GO:0000003",
            "go_id": "GO:0000003",
            "term_name": "leaf function",
        },
    ]


def test_go_db_select_terms_keeps_one_row_per_term_for_subset(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    df_terms = GODatabase.from_obo(file_in).select_terms(subset_id="goslim_generic")

    assert df_terms["go_id"].to_list() == [
        "GO:0000001",
        "GO:0000002",
        "GO:0005575",
        "GO:0005737",
    ]
    assert df_terms["subset_id"].to_list() == ["goslim_generic"] * 4


def test_go_db_extracts_subcell_from_cellular_component(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    df_subcell = GODatabase.from_obo(file_in).extract_subcell()

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

    df_subcell_with_obsolete = GODatabase.from_obo(file_in).extract_subcell(
        include_obsolete=True
    )
    assert df_subcell_with_obsolete["go_id"].to_list() == [
        "GO:0000004",
        "GO:0005575",
        "GO:0005737",
    ]


def test_go_db_writes_subcell_parquet(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    path = tmp_path / "subcell.parquet"
    write_minimal_obo(file_in)

    path_written = GODatabase.from_obo(file_in).write_subcell(path)

    assert path_written == path
    assert pl.read_parquet(path).height == 2
