import gzip
import tarfile
import zipfile
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract import GODatabase
from bioextract.go.go import GoSubsetId


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
        df_term = pl.read_database(  # pyright: ignore[reportUnknownMemberType]  # Polars-DuckDB boundary
            "SELECT * FROM term", connection
        )
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


def test_go_duckdb_reopen_preserves_domain_and_native_sql_behavior(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    source = GODatabase.from_obo(source_path)
    expected_terms = source.select_terms(
        term_ids=["GO:1234567", "GO:0005575"],
        include_obsolete=True,
    )
    expected_subsets = source.list_subsets()
    source.write_duckdb(publication)

    reopened = GODatabase.from_duckdb(publication)
    assert reopened.select_terms(
        term_ids=["GO:1234567", "GO:0005575"],
        include_obsolete=True,
    ).equals(expected_terms)
    assert reopened.list_subsets().equals(expected_subsets)

    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        assert first.execute(
            "SELECT term_name FROM term WHERE go_id='GO:0005737'"
        ).fetchone() == ("cytoplasm",)
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_go_duckdb_reopen_uses_declared_schemas_for_nullable_synonym_types(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    synonym_lines = [
        f'synonym: "untyped synonym {index:03d}" EXACT []' for index in range(100)
    ]
    synonym_lines.append('synonym: "typed synonym" EXACT systematic_synonym []')
    with source_path.open("a", encoding="utf-8") as source_file:
        source_file.write(
            "\n[Term]\n"
            "id: GO:0000005\n"
            "name: synonym inference fixture\n"
            "namespace: biological_process\n" + "\n".join(synonym_lines) + "\n"
        )

    source = GODatabase.from_obo(source_path)
    expected_tidy = {
        name: frame.collect() for name, frame in source.build_tidy().frames.items()
    }
    expected_terms = source.select_terms(term_ids=["GO:0000005"])
    expected_subsets = source.list_subsets()
    source.write_duckdb(publication)

    reopened = GODatabase.from_duckdb(publication)
    actual_tidy = reopened.build_tidy()
    assert set(actual_tidy.frames) == set(expected_tidy)
    for name, expected_frame in expected_tidy.items():
        assert actual_tidy.frames[name].collect().equals(expected_frame)
    assert reopened.select_terms(term_ids=["GO:0000005"]).equals(expected_terms)
    assert reopened.list_subsets().equals(expected_subsets)
    assert (
        actual_tidy.frames["synonym"]
        .filter(pl.col("synonym_type_name") == "systematic_synonym")
        .select("synonym_scope")
        .collect()
        .item()
        == "EXACT"
    )


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


def test_go_db_selects_ancestors_and_reports_unmatched(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    selection = GODatabase.from_obo(file_in).select_ancestors(
        [" GO:1234567 ", "GO:0000003", "GO:0000004", "GO:9999999"],
        target_subset_id=GoSubsetId.GOSLIM_GENERIC,
        include_self=True,
    )

    assert selection.extract_ancestors().to_dicts() == [
        {
            "input_go_id": "GO:1234567",
            "go_id": "GO:0000002",
            "ancestor_go_id": "GO:0000002",
            "ancestor_term_name": "child process",
            "ancestor_namespace": "biological_process",
            "min_distance": 0,
            "target_subset_id": "goslim_generic",
        },
        {
            "input_go_id": "GO:1234567",
            "go_id": "GO:0000002",
            "ancestor_go_id": "GO:0000001",
            "ancestor_term_name": "root process",
            "ancestor_namespace": "biological_process",
            "min_distance": 1,
            "target_subset_id": "goslim_generic",
        },
        {
            "input_go_id": "GO:0000003",
            "go_id": "GO:0000003",
            "ancestor_go_id": "GO:0000002",
            "ancestor_term_name": "child process",
            "ancestor_namespace": "biological_process",
            "min_distance": 1,
            "target_subset_id": "goslim_generic",
        },
        {
            "input_go_id": "GO:0000003",
            "go_id": "GO:0000003",
            "ancestor_go_id": "GO:0000001",
            "ancestor_term_name": "root process",
            "ancestor_namespace": "biological_process",
            "min_distance": 2,
            "target_subset_id": "goslim_generic",
        },
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {
            "input_go_id": "GO:0000004",
            "resolved_go_id": "GO:0000004",
            "reason": "obsolete_excluded",
        },
        {
            "input_go_id": "GO:9999999",
            "resolved_go_id": None,
            "reason": "not_found",
        },
    ]


def test_go_db_ancestor_selection_handles_no_match_and_invalid_input(
    tmp_path: Path,
) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)
    database = GODatabase.from_obo(file_in)

    unmatched = database.select_ancestors(
        ["GO:0000001"],
        target_subset_id="missing_subset",
    ).extract_unmatched_ids()
    assert unmatched.to_dicts() == [
        {
            "input_go_id": "GO:0000001",
            "resolved_go_id": "GO:0000001",
            "reason": "no_matching_ancestor",
        }
    ]

    with pytest.raises(ValueError, match="Invalid GO identifier"):
        database.select_ancestors(["not-a-go-id"])


def test_go_ancestor_selection_empty_schema(tmp_path: Path) -> None:
    file_in = tmp_path / "go-basic.obo"
    write_minimal_obo(file_in)

    selection = GODatabase.from_obo(file_in).select_ancestors([])

    assert selection.extract_ancestors().schema == {
        "input_go_id": pl.String,
        "go_id": pl.String,
        "ancestor_go_id": pl.String,
        "ancestor_term_name": pl.String,
        "ancestor_namespace": pl.String,
        "min_distance": pl.Int32,
        "target_subset_id": pl.String,
    }
    assert selection.extract_unmatched_ids().schema == {
        "input_go_id": pl.String,
        "resolved_go_id": pl.String,
        "reason": pl.Enum(["not_found", "obsolete_excluded", "no_matching_ancestor"]),
    }


@pytest.mark.parametrize(
    (
        "term_ids",
        "target_subset_id",
        "include_self",
        "resolve_alt_ids",
        "include_obsolete",
    ),
    [
        (["GO:1234567"], None, False, True, False),
        (["GO:1234567"], "goslim_generic", True, False, False),
        (["GO:0000004"], "goslim_generic", True, True, True),
        ([], "goslim_generic", True, True, False),
    ],
)
def test_go_ancestor_selection_policy_parity(
    tmp_path: Path,
    term_ids: list[str],
    target_subset_id: str | None,
    include_self: bool,
    resolve_alt_ids: bool,
    include_obsolete: bool,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    source = GODatabase.from_obo(source_path)
    source.write_duckdb(publication)
    published = GODatabase.from_duckdb(publication)

    source_selection = source.select_ancestors(
        term_ids,
        target_subset_id=target_subset_id,
        include_self=include_self,
        resolve_alt_ids=resolve_alt_ids,
        include_obsolete=include_obsolete,
    )
    publication_selection = published.select_ancestors(
        term_ids,
        target_subset_id=target_subset_id,
        include_self=include_self,
        resolve_alt_ids=resolve_alt_ids,
        include_obsolete=include_obsolete,
    )
    assert publication_selection.extract_ancestors().equals(
        source_selection.extract_ancestors()
    )
    assert publication_selection.extract_unmatched_ids().equals(
        source_selection.extract_unmatched_ids()
    )


def test_go_duckdb_ancestor_selection_matches_source_and_avoids_full_frame_loader(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    source = GODatabase.from_obo(source_path)
    source.write_duckdb(publication)
    published = GODatabase.from_duckdb(publication)

    def fail_class_full_frame_loader(
        database: GODatabase, *args: object, **kwargs: object
    ) -> None:
        del database, args, kwargs
        raise AssertionError("domain query loaded complete GO frames")

    monkeypatch.setattr(
        GODatabase,
        "_read_publication_frames",
        fail_class_full_frame_loader,
    )

    source_selection = source.select_ancestors(
        ["GO:1234567", "GO:0000003", "GO:9999999"],
        target_subset_id="goslim_generic",
        include_self=True,
    )
    publication_selection = published.select_ancestors(
        ["GO:1234567", "GO:0000003", "GO:9999999"],
        target_subset_id="goslim_generic",
        include_self=True,
    )
    assert publication_selection.extract_ancestors().equals(
        source_selection.extract_ancestors()
    )
    assert publication_selection.extract_unmatched_ids().equals(
        source_selection.extract_unmatched_ids()
    )


def test_go_publication_domain_queries_do_not_use_full_frame_loader(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    source = GODatabase.from_obo(source_path)
    source.write_duckdb(publication)
    published = GODatabase.from_duckdb(publication)

    def fail_class_full_frame_loader(
        database: GODatabase, *args: object, **kwargs: object
    ) -> None:
        del database, args, kwargs
        raise AssertionError("domain query loaded complete GO frames")

    monkeypatch.setattr(
        GODatabase,
        "_read_publication_frames",
        fail_class_full_frame_loader,
    )
    assert published.select_terms(term_ids=["GO:1234567"]).height == 1
    assert published.list_subsets().height == 2


def test_go_ancestor_selection_reuses_one_publication_resolution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_path = tmp_path / "go-basic.obo"
    publication = tmp_path / "go.duckdb"
    write_minimal_obo(source_path)
    source = GODatabase.from_obo(source_path)
    source.write_duckdb(publication)
    published = GODatabase.from_duckdb(publication)

    calls = 0
    connect = GODatabase.connect

    def counted_connect(database: GODatabase) -> duckdb.DuckDBPyConnection:
        nonlocal calls
        calls += 1
        return connect(database)

    monkeypatch.setattr(GODatabase, "connect", counted_connect)
    selection = published.select_ancestors(["GO:1234567"])
    selection.extract_ancestors()
    selection.extract_unmatched_ids()
    assert calls == 1
