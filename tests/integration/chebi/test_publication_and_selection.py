from __future__ import annotations

import gzip
import tarfile
import zipfile
from pathlib import Path

import duckdb
import pytest

import bioextract.chebi._query as chebi_query
from bioextract import ChEBIDatabase

COMPOUNDS = """id\tname\tstatus_id\tsource\tparent_id\tmerge_type\tchebi_accession\tdefinition\tascii_name\tstars\tmodified_on\trelease_date
1\twater\t1\tChEBI\t\t\tCHEBI:1\twater definition\twater\t3\t2026-01-01\t2026-01-01
"""
NAMES = """id\tcompound_id\tname\ttype\tstatus_id\tadapted\tlanguage_code\tascii_name
10\t1\taqua\tSYNONYM\t1\tF\ten\taqua
"""
OBO = """format-version: 1.2

[Term]
id: CHEBI:1
name: water
def: "water definition" []
synonym: "aqua" EXACT []
xref: KEGG:C00001
is_a: CHEBI:0 ! chemical entity
relationship: has_role CHEBI:2
"""


def _write_release(directory: Path, *, gzip_compounds: bool = False) -> None:
    directory.mkdir()
    if gzip_compounds:
        with gzip.open(
            directory / "compounds.tsv.gz",
            "wt",
            encoding="utf-8",
        ) as handle:
            handle.write(COMPOUNDS)
    else:
        (directory / "compounds.tsv").write_text(COMPOUNDS, encoding="utf-8")
    (directory / "names.tsv").write_text(NAMES, encoding="utf-8")


@pytest.mark.parametrize("archive_kind", ["zip", "tar"])
def test_release_archive_is_detected_from_content(
    tmp_path: Path,
    archive_kind: str,
) -> None:
    file_archive = tmp_path / "release.snapshot"
    if archive_kind == "zip":
        with zipfile.ZipFile(file_archive, "w") as archive:
            archive.writestr("nested/compounds.tsv", COMPOUNDS)
            archive.writestr("nested/names.tsv", NAMES)
    else:
        directory = tmp_path / "release"
        _write_release(directory)
        with tarfile.open(file_archive, "w") as archive:
            archive.add(directory / "compounds.tsv", arcname="nested/compounds.tsv")
            archive.add(directory / "names.tsv", arcname="nested/names.tsv")

    result = ChEBIDatabase.from_release(file_archive).write_duckdb(
        tmp_path / f"{archive_kind}.duckdb"
    )
    assert result.row_counts == {"compound": 1, "compound_name": 1}


@pytest.mark.parametrize("container", ["plain", "gzip", "zip", "tar"])
def test_obo_input_container_is_detected_from_content(
    tmp_path: Path,
    container: str,
) -> None:
    file_obo = tmp_path / f"chebi-{container}.snapshot"
    if container == "plain":
        file_obo.write_text(OBO, encoding="utf-8")
    elif container == "gzip":
        with gzip.open(file_obo, "wt", encoding="utf-8") as handle:
            handle.write(OBO)
    elif container == "zip":
        with zipfile.ZipFile(file_obo, "w") as archive:
            archive.writestr("ontology/chebi.obo", OBO)
    else:
        file_plain = tmp_path / "chebi.obo"
        file_plain.write_text(OBO, encoding="utf-8")
        with tarfile.open(file_obo, "w") as archive:
            archive.add(file_plain, arcname="ontology/chebi.obo")

    path = tmp_path / f"{container}.duckdb"
    result = ChEBIDatabase.from_obo(file_obo).write_duckdb(path)
    assert result.row_counts["compound"] == 1
    assert result.row_counts["compound_relation"] == 0
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT name, scope FROM compound_name"
        ).fetchone() == ("aqua", "EXACT")
        assert connection.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone() == (2,)


def test_chemont_relations_share_container_but_not_domain_tables(
    tmp_path: Path,
) -> None:
    file_chebi = tmp_path / "chebi.obo"
    file_chemont = tmp_path / "chemont.obo"
    file_chebi.write_text(OBO, encoding="utf-8")
    file_chemont.write_text(
        OBO.replace("CHEBI:", "CHEMONT:"),
        encoding="utf-8",
    )

    result = ChEBIDatabase.from_obo(
        file_chebi,
        chemont_obo=file_chemont,
    ).write_duckdb(tmp_path / "combined.duckdb")

    assert "compound" in result.tables
    assert "chemont_term" in result.tables
    assert "chebi_chemont" not in result.tables


DOMAIN_OBO = """format-version: 1.2

[Term]
id: CHEBI:1
name: water
def: "water definition" []
alt_id: CHEBI:100
subset: 3_STAR
synonym: "aqua" EXACT []
xref: kegg.compound:C00001
xref: HMDB:HMDB0002111
property_value: http://purl.obolibrary.org/obo/chebi/formula "H2O" xsd:string
property_value: http://purl.obolibrary.org/obo/chebi/inchi "InChI=1S/H2O/h1H2" xsd:string
property_value: http://purl.obolibrary.org/obo/chebi/inchikey "XLYOFNOQVPJJNP-UHFFFAOYSA-N" xsd:string
is_a: CHEBI:2 ! parent

[Term]
id: CHEBI:2
name: parent
subset: 2_STAR

[Term]
id: CHEBI:3
name: obsolete
subset: 3_STAR
is_obsolete: true

[Term]
id: CHEBI:4
name: low star
subset: 1_STAR

[Term]
id: CHEBI:5
name: orphan source
subset: 3_STAR
alt_id: CHEBI:500
is_a: CHEBI:999 ! missing target
"""


def _write_domain_publication(tmp_path: Path) -> Path:
    file_obo = tmp_path / "chebi.obo"
    file_obo.write_text(DOMAIN_OBO, encoding="utf-8")
    file_sdf = tmp_path / "chebi.sdf"
    file_sdf.write_text(
        """water
  bioextract

  0  0  0  0  0  0            999 V2000
M  END
> <ChEBI ID>
CHEBI:1

$$$$
orphan
  bioextract

  0  0  0  0  0  0            999 V2000
M  END
> <ChEBI ID>
CHEBI:999

$$$$
""",
        encoding="utf-8",
    )
    path = tmp_path / "chebi.duckdb"
    ChEBIDatabase.from_obo(file_obo, sdf=file_sdf).write_duckdb(path)
    return path


def test_domain_selection_supports_shared_ids_namespaces_and_relations(
    tmp_path: Path,
) -> None:
    path = _write_domain_publication(tmp_path)
    database = ChEBIDatabase.from_duckdb(path)

    assert database.select_compounds(
        ["1", "CHEBI:100"], namespace="chebi"
    ).extract_matches()["chebi_id"].to_list() == ["CHEBI:1", "CHEBI:1"]
    assert database.select_compounds(
        ["InChI=1S/H2O/h1H2"], namespace="inchi"
    ).extract_compounds().select("chebi_id", "preferred_name").to_dicts() == [
        {"chebi_id": "CHEBI:1", "preferred_name": "water"}
    ]
    assert database.select_compounds(
        ["C00001"], namespace="kegg.compound"
    ).extract_matches()["chebi_id"].to_list() == ["CHEBI:1"]
    selection = database.select_groups(
        {"first": ["CHEBI:1"], "second": ["CHEBI:2"]},
        namespace="chebi",
    )
    assert selection.extract_ancestors().select(
        "group_id", "chebi_id", "ancestor_chebi_id"
    ).to_dicts() == [
        {
            "group_id": "first",
            "chebi_id": "CHEBI:1",
            "ancestor_chebi_id": "CHEBI:2",
        }
    ]
    assert (
        database.select_compounds(["CHEBI:1"], namespace="chebi")
        .extract_structures()["molfile"]
        .str.contains("M  END")
        .all()
    )

    with database.connect() as connection:
        assert connection.execute(
            "SELECT preferred_name FROM compound WHERE chebi_id='CHEBI:1'"
        ).fetchone() == ("water",)
        with pytest.raises(duckdb.Error):
            connection.execute("CREATE TABLE forbidden(value INTEGER)")


def test_grouped_selection_resolves_unique_ids_once_and_reuses_candidates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    database = ChEBIDatabase.from_duckdb(_write_domain_publication(tmp_path))
    input_table_calls: list[tuple[tuple[str, str], ...]] = []
    original_create_input_table = chebi_query._create_input_table  # pyright: ignore[reportPrivateUsage]

    def counted_create_input_table(
        connection: duckdb.DuckDBPyConnection,
        rows: tuple[chebi_query._InputRow, ...],  # pyright: ignore[reportPrivateUsage]
    ) -> None:
        input_table_calls.append(
            tuple((row.input_id, row.lookup_value) for row in rows)
        )
        original_create_input_table(connection, rows)

    monkeypatch.setattr(
        chebi_query,
        "_create_input_table",
        counted_create_input_table,
    )
    selection = database.select_groups(
        {
            " first ": ["1", "CHEBI:1", "CHEBI:404"],
            "second": ["CHEBI:1", "404"],
            "empty": [],
        },
        namespace="chebi",
    )

    assert selection._group_ids == (  # pyright: ignore[reportPrivateUsage]
        "empty",
        "first",
        "second",
    )
    assert selection.extract_matches().select(
        "group_id", "input_id", "chebi_id"
    ).to_dicts() == [
        {"group_id": "first", "input_id": "CHEBI:1", "chebi_id": "CHEBI:1"},
        {"group_id": "second", "input_id": "CHEBI:1", "chebi_id": "CHEBI:1"},
    ]
    selection.extract_compounds()
    selection.extract_matches()
    assert selection.extract_unmatched_ids().select(
        "group_id", "input_id", "reason"
    ).to_dicts() == [
        {"group_id": "first", "input_id": "CHEBI:404", "reason": "not_found"},
        {"group_id": "second", "input_id": "CHEBI:404", "reason": "not_found"},
    ]
    selection.extract_unmatched_ids()
    assert input_table_calls == [
        (
            ("CHEBI:1", "CHEBI:1"),
            ("CHEBI:404", "CHEBI:404"),
        ),
        (
            ("CHEBI:1", "CHEBI:1"),
            ("CHEBI:404", "CHEBI:404"),
        ),
    ]


def test_unmatched_reason_precedence_and_dynamic_namespace_error(
    tmp_path: Path,
) -> None:
    database = ChEBIDatabase.from_duckdb(_write_domain_publication(tmp_path))
    reasons = (
        database.select_compounds(
            ["CHEBI:3", "CHEBI:4", "CHEBI:999", "CHEBI:404"],
            namespace="chebi",
            min_star_rating=2,
        )
        .extract_unmatched_ids()
        .sort("input_id")
        .to_dicts()
    )
    assert reasons == [
        {
            "input_id": "CHEBI:3",
            "input_namespace": "chebi",
            "reason": "obsolete_excluded",
        },
        {
            "input_id": "CHEBI:4",
            "input_namespace": "chebi",
            "reason": "below_min_star_rating",
        },
        {
            "input_id": "CHEBI:404",
            "input_namespace": "chebi",
            "reason": "not_found",
        },
        {
            "input_id": "CHEBI:999",
            "input_namespace": "chebi",
            "reason": "invalid_canonical_target",
        },
    ]
    with pytest.raises(ValueError, match=r"available:.*kegg\.compound"):
        database.select_compounds(["x"], namespace="unknown")
