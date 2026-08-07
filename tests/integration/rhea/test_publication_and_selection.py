from __future__ import annotations

import gzip
import zipfile
from pathlib import Path

import duckdb
import polars as pl
import pytest

import bioextract.rhea._query as rhea_query
from bioextract import RheaDatabase
from bioextract.rhea.constant import RheaNamespace

RDF_FIXTURE = """\
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
  xmlns:rh="http://rdf.rhea-db.org/">
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000">
    <rh:accession>RHEA:10000</rh:accession>
    <rh:equation>A = B</rh:equation>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
    <rh:isChemicallyBalanced>true</rh:isChemicallyBalanced>
    <rh:isTransport>false</rh:isTransport>
    <rh:side rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:side rdf:resource="http://rdf.rhea-db.org/10000_R"/>
    <rh:citation rdf:resource="http://rdf.rhea-db.org/12345"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10001">
    <rh:accession>RHEA:10001</rh:accession>
    <rh:substrates rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:products rdf:resource="http://rdf.rhea-db.org/10000_R"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10002">
    <rh:accession>RHEA:10002</rh:accession>
    <rh:substrates rdf:resource="http://rdf.rhea-db.org/10000_R"/>
    <rh:products rdf:resource="http://rdf.rhea-db.org/10000_L"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10003">
    <rh:accession>RHEA:10003</rh:accession>
    <rh:substratesOrProducts rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:substratesOrProducts rdf:resource="http://rdf.rhea-db.org/10000_R"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10004">
    <rh:accession>RHEA:10004</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Obsolete"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_L">
    <rh:curatedOrder>1</rh:curatedOrder>
    <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_A"/>
    <rh:contains1 rdf:resource="http://rdf.rhea-db.org/Participant_A"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_R">
    <rh:curatedOrder>2</rh:curatedOrder>
    <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_B"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/contains1">
    <rh:coefficient>2</rh:coefficient>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_A">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_A"/>
    <rh:location rdf:resource="http://rdf.rhea-db.org/Cytoplasm"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_B">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_B"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_A">
    <rh:id>1</rh:id>
    <rh:accession>CHEBI:1</rh:accession>
    <rh:name>A</rh:name>
    <rh:formula>C</rh:formula>
    <rh:charge>-1</rh:charge>
    <rh:chebi rdf:resource="http://purl.obolibrary.org/obo/CHEBI_1"/>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/SmallMolecule"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_B">
    <rh:id>2</rh:id>
    <rh:accession>GENERIC:2</rh:accession>
    <rh:name>B</rh:name>
    <rh:charge>(-1)n</rh:charge>
    <rh:reactivePart rdf:resource="http://rdf.rhea-db.org/Compound_Part"/>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/Polymer"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_Part">
    <rh:position>C1</rh:position>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/ReactivePart"/>
  </rdf:Description>
</rdf:RDF>
"""

SDF_FIXTURE = """\
CHEBI:1
  bioextract

  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <ACCESSION>
CHEBI:1

> <ROLE>
small molecule

> <CHEBI_XREF>
CHEBI:1

> <Formula>
C

> <Charge>
(-1)n

> <Rhea_ascii_name>
A

$$$$
"""


def _write(path: Path, value: str) -> Path:
    path.write_text(value, encoding="utf-8")
    return path


def _write_release(root: Path) -> Path:
    raw = root / "raw"
    raw.mkdir(parents=True)
    _write(raw / "LICENSE.txt", "Creative Commons\n")
    _write(raw / "chebiId_name.tsv", "CHEBI:1\t A\n")
    _write(
        raw / "chebi_pH7_3_mapping.tsv",
        "CHEBI\tCHEBI_PH7_3\tORIGIN\n1\t1\tcomputation\n",
    )
    _write(
        raw / "rhea-directions.tsv",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "10000\t10001\t10002\t10003\n",
    )
    _write(raw / "rhea-obsoletes.tsv", "RHEA_ID\n")
    _write(raw / "rhea-reaction-smiles.tsv", "10000\tA>>B\n")
    _write(
        raw / "rhea-relationships.tsv",
        "FROM_REACTION_ID\tTYPE\tTO_REACTION_ID\n10000\tis_a\t10004\n",
    )
    _write(
        raw / "rhea-release.properties",
        "rhea.release.number=141\nrhea.release.date=2026-06-10\n",
    )
    _write(raw / "rhea.rdf", RDF_FIXTURE)
    _write(raw / "rhea.sdf", SDF_FIXTURE)
    header = "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n"
    _write(raw / "rhea2ec.tsv", header + "10000\tUN\t10000\t1.1.1.1\n")
    _write(raw / "rhea2go.tsv", header + "10000\tUN\t10000\tGO:0000001\n")
    _write(
        raw / "rhea2uniprot_sprot.tsv",
        header + "10000\tUN\t10000\tP00001\n",
    )
    _write(
        raw / "rhea2uniprot_trembl.tsv",
        header + "10001\tLR\t10000\tA00001\n",
    )
    _write(
        raw / "rhea2xrefs.tsv",
        "RHEA_ID\tDIRECTION\tMASTER_ID\tID\tDB\n"
        "10000\tUN\t10000\t1.1.1.1\tEC\n"
        "10000\tUN\t10000\tGO:0000001\tGO\n"
        "10000\tUN\t10000\tRXN-1\tECOCYC\n"
        "10000\tUN\t10000\tR00001\tKEGG_REACTION\n"
        "10000\tUN\t10000\tM0001\tMACIE\n"
        "10000\tUN\t10000\tMETA:RXN-1\tMETACYC\n"
        "10000\tUN\t10000\tR-HSA-1\tREACTOME\n",
    )
    return root


def test_reaction_files_write_semantic_tables(tmp_path: Path) -> None:
    file_rdf = _write(tmp_path / "rhea.rdf", RDF_FIXTURE)
    file_directions = _write(
        tmp_path / "directions.tsv",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "10000\t10001\t10002\t10003\n",
    )
    path = tmp_path / "rhea.duckdb"

    report = RheaDatabase.from_files(
        rdf=file_rdf,
        directions=file_directions,
    ).write_duckdb(path)

    assert report.scope == "reactions"
    assert report.row_counts["reaction"] == 5
    assert "compound_structure" not in report.tables
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT master_id, direction FROM reaction WHERE rhea_id = 10002"
        ).fetchone() == (10000, "RL")
        assert connection.execute(
            """
            SELECT coefficient_text, coefficient_numeric, location
            FROM reaction_participant
            WHERE participant_id = 'Participant_A'
            """
        ).fetchone() == ("2", 2.0, "cytoplasm")
        assert connection.execute(
            """
            SELECT directional_role FROM reaction_participant_direction
            WHERE rhea_id = 10002 AND participant_id = 'Participant_A'
            """
        ).fetchone() == ("product",)
        assert connection.execute(
            "SELECT charge_text, charge_numeric FROM compound "
            "WHERE compound_id = 'Compound_B'"
        ).fetchone() == ("(-1)n", None)
        assert connection.execute(
            "SELECT master_id, direction, is_obsolete "
            "FROM reaction WHERE rhea_id = 10004"
        ).fetchone() == (None, None, True)


def test_gzip_is_detected_from_content_not_suffix(tmp_path: Path) -> None:
    file_rdf = tmp_path / "rhea.data"
    with gzip.open(file_rdf, "wt", encoding="utf-8") as handle:
        handle.write(RDF_FIXTURE)
    file_directions = _write(
        tmp_path / "directions.data",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "10000\t10001\t10002\t10003\n",
    )

    report = RheaDatabase.from_files(
        rdf=file_rdf,
        directions=file_directions,
    ).write_duckdb(tmp_path / "gzip.duckdb")

    assert report.row_counts["reaction"] == 5


def test_participant_constructor_creates_only_supplied_tables(
    tmp_path: Path,
) -> None:
    file_names = tmp_path / "names.data"
    with gzip.open(file_names, "wt", encoding="utf-8") as handle:
        handle.write("CHEBI:1\t A\nCHEBI:2\t B\n")
    path = tmp_path / "participants.duckdb"

    report = RheaDatabase.from_files(chebi_names=file_names).write_duckdb(path)

    assert report.tables == ("chebi_name",)
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT name FROM chebi_name ORDER BY chebi_id"
        ).fetchall() == [("A",), ("B",)]


def test_cross_reference_constructor_builds_views(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    raw = release / "raw"
    path = tmp_path / "xrefs.duckdb"

    RheaDatabase.from_files(
        xrefs=raw / "rhea2xrefs.tsv",
        uniprot_sprot=raw / "rhea2uniprot_sprot.tsv",
    ).write_duckdb(path)

    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute("SELECT ec_number FROM reaction_ec").fetchone() == (
            "1.1.1.1",
        )
        assert connection.execute(
            "SELECT uniprot_section FROM reaction_uniprot"
        ).fetchone() == ("Swiss-Prot",)


def test_complete_release_directory_and_archive(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    file_directory = tmp_path / "directory.duckdb"
    directory_report = RheaDatabase.from_files(release).write_duckdb(file_directory)

    assert directory_report.release_number == 141
    assert directory_report.release_date == "2026-06-10"
    assert directory_report.row_counts["compound_structure"] == 1
    assert directory_report.row_counts["reaction_uniprot"] == 2
    with duckdb.connect(str(file_directory), read_only=True) as connection:
        assert connection.execute(
            "SELECT charge_text, charge_numeric FROM compound_structure"
        ).fetchone() == ("(-1)n", None)

    file_archive = tmp_path / "release.zip"
    with zipfile.ZipFile(file_archive, mode="w") as archive:
        for path in release.rglob("*"):
            if path.is_file():
                archive.write(path, path.relative_to(tmp_path))
    file_archived_db = tmp_path / "archive.duckdb"
    archive_report = RheaDatabase.from_files(file_archive).write_duckdb(
        file_archived_db
    )

    assert archive_report.row_counts == directory_report.row_counts
    assert all(
        str(file_archive) in value for value in archive_report.source_files.values()
    )


def test_from_duckdb_selects_exact_reactions_and_domain_relations(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)

    database = RheaDatabase.from_duckdb(path)
    selection = database.select_reactions(["CHEBI:1"], namespace="chebi")

    assert selection.extract_matches().select(
        "input_id", "input_namespace", "rhea_id", "direction"
    ).to_dicts() == [
        {
            "input_id": "CHEBI:1",
            "input_namespace": "chebi",
            "rhea_id": 10000,
            "direction": "UN",
        },
        {
            "input_id": "CHEBI:1",
            "input_namespace": "chebi",
            "rhea_id": 10001,
            "direction": "LR",
        },
        {
            "input_id": "CHEBI:1",
            "input_namespace": "chebi",
            "rhea_id": 10002,
            "direction": "RL",
        },
        {
            "input_id": "CHEBI:1",
            "input_namespace": "chebi",
            "rhea_id": 10003,
            "direction": "BI",
        },
    ]
    participant_roles = (
        selection.extract_participants()
        .filter(pl.col("participant_id") == "Participant_A")
        .select("rhea_id", "side", "directional_role", "chebi_id")
        .to_dicts()
    )
    assert participant_roles == [
        {
            "rhea_id": 10000,
            "side": "L",
            "directional_role": None,
            "chebi_id": "CHEBI:1",
        },
        {
            "rhea_id": 10001,
            "side": "L",
            "directional_role": "substrate",
            "chebi_id": "CHEBI:1",
        },
        {
            "rhea_id": 10002,
            "side": "L",
            "directional_role": "product",
            "chebi_id": "CHEBI:1",
        },
        {
            "rhea_id": 10003,
            "side": "L",
            "directional_role": None,
            "chebi_id": "CHEBI:1",
        },
    ]
    assert selection.extract_reactions().filter(pl.col("rhea_id") == 10000)[
        "reaction_smiles"
    ].to_list() == ["A>>B"]
    assert selection.extract_publications()["pubmed_id"].to_list() == ["12345"]
    assert selection.extract_relationships()["relation_type"].unique().to_list() == [
        "is_a"
    ]


def test_grouped_selection_preserves_lineage_and_unmatched_ids(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)

    selection = RheaDatabase.from_duckdb(path).select_groups(
        {
            "known": ["1.1.1.1"],
            "mixed": ["1.1.1.1", "9.9.9.9"],
        },
        namespace="ec",
    )

    assert selection.extract_matches().select(
        "group_id", "input_id", "rhea_id"
    ).to_dicts() == [
        {"group_id": "known", "input_id": "1.1.1.1", "rhea_id": 10000},
        {"group_id": "mixed", "input_id": "1.1.1.1", "rhea_id": 10000},
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {
            "group_id": "mixed",
            "input_id": "9.9.9.9",
            "input_namespace": "ec",
        }
    ]


def test_grouped_selection_resolves_unique_ids_once_and_reuses_mapping(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)
    input_table_calls: list[tuple[tuple[str, str], ...]] = []
    original_create_input_table = rhea_query._create_input_table  # pyright: ignore[reportPrivateUsage]

    def counted_create_input_table(
        connection: duckdb.DuckDBPyConnection,
        rows: tuple[rhea_query._InputRow, ...],  # pyright: ignore[reportPrivateUsage]
    ) -> None:
        input_table_calls.append(
            tuple((row.input_id, row.lookup_value) for row in rows)
        )
        original_create_input_table(connection, rows)

    monkeypatch.setattr(
        rhea_query,
        "_create_input_table",
        counted_create_input_table,
    )
    selection = RheaDatabase.from_duckdb(path).select_groups(
        {
            " case ": ["10000", "RHEA:10000", "99999"],
            "control": ["RHEA:10000", "99999"],
            "empty": [],
        },
        namespace="rhea",
    )

    assert selection._group_ids == (  # pyright: ignore[reportPrivateUsage]
        "case",
        "control",
        "empty",
    )
    assert selection.extract_matches().select(
        "group_id", "input_id", "rhea_id"
    ).to_dicts() == [
        {"group_id": "case", "input_id": "RHEA:10000", "rhea_id": 10000},
        {"group_id": "control", "input_id": "RHEA:10000", "rhea_id": 10000},
    ]
    selection.extract_reactions()
    selection.extract_matches()
    assert selection.extract_unmatched_ids().select(
        "group_id", "input_id"
    ).to_dicts() == [
        {"group_id": "case", "input_id": "RHEA:99999"},
        {"group_id": "control", "input_id": "RHEA:99999"},
    ]
    assert input_table_calls == [
        (
            ("RHEA:10000", "10000"),
            ("RHEA:99999", "99999"),
        )
    ]


@pytest.mark.parametrize(
    ("namespace", "input_id"),
    [
        ("rhea", "RHEA:10000"),
        ("chebi", "CHEBI:1"),
        ("uniprot", "sp|P00001|TEST_HUMAN"),
        ("ec", "1.1.1.1"),
        ("go", "go:0000001"),
        ("ecocyc", "RXN-1"),
        ("kegg_reaction", "R00001"),
        ("macie", "M0001"),
        ("metacyc", "META:RXN-1"),
        ("reactome", "R-HSA-1"),
    ],
)
def test_supported_namespaces_resolve_official_mappings(
    tmp_path: Path,
    namespace: RheaNamespace,
    input_id: str,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)

    matches = (
        RheaDatabase.from_duckdb(path)
        .select_reactions(
            [input_id],
            namespace=namespace,
        )
        .extract_matches()
    )

    assert 10000 in matches["rhea_id"].to_list()


def test_obsolete_policy_is_explicit(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)
    database = RheaDatabase.from_duckdb(path)

    excluded = database.select_reactions(["10004"], namespace="rhea")
    included = database.select_reactions(
        ["RHEA:10004"],
        namespace="rhea",
        include_obsolete=True,
    )

    assert excluded.extract_matches().is_empty()
    assert excluded.extract_unmatched_ids()["input_id"].to_list() == ["RHEA:10004"]
    assert included.extract_matches()["rhea_id"].to_list() == [10004]
