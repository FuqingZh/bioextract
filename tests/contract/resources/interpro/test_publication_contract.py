from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import pytest

from bioextract.errors import CapabilityError, IntegrityError
from bioextract.interpro import InterProDatabase


def _source(
    tmp_path: Path,
    *,
    with_xml: bool = True,
    mapping_rows: tuple[str, ...] = ("P12345\tIPR000001\tKringle\tPF00051\t10\t80",),
) -> InterProDatabase:
    mapping = tmp_path / "protein2ipr.dat.gz"
    with gzip.open(mapping, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(mapping_rows) + "\n")
    if not with_xml:
        return InterProDatabase.from_mapping_files(protein_to_interpro=mapping)
    xml = tmp_path / "interpro.xml.gz"
    with gzip.open(xml, "wt", encoding="utf-8") as handle:
        handle.write(
            """<interprodb>
<release><dbinfo dbname="INTERPRO" version="108.0"/></release>
<interpro id="IPR000001" type="Domain">
<name>Kringle</name><member_list>
<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>
</member_list></interpro></interprodb>"""
        )
    return InterProDatabase.from_mapping_files(
        protein_to_interpro=mapping,
        interpro_xml=xml,
    )


@pytest.mark.parametrize(
    ("with_xml", "profile", "capabilities", "tables", "source_roles"),
    [
        (
            False,
            "interpro-protein2ipr-v1",
            "mapping",
            {"mapping": "canonical"},
            {"protein_to_interpro"},
        ),
        (
            True,
            "interpro-protein2ipr-xml-v1",
            "mapping,pfam",
            {
                "mapping": "canonical",
                "protein_term": "compact",
                "term": "compact",
                "term_xref": "compact",
            },
            {"protein_to_interpro", "interpro_xml"},
        ),
    ],
)
def test_capability_publication_has_exact_metadata_v1_inventory_and_reopens(
    tmp_path: Path,
    with_xml: bool,
    profile: str,
    capabilities: str,
    tables: dict[str, str],
    source_roles: set[str],
) -> None:
    path = tmp_path / "interpro.duckdb"
    result = _source(tmp_path, with_xml=with_xml).write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)
    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        metadata = dict(
            first.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "1"
        assert metadata["bioextract.resource_name"] == "interpro"
        assert metadata["bioextract.resource_schema_version"] == "interpro-v1"
        assert metadata["bioextract.source_schema_profile"] == profile
        assert metadata["bioextract.capabilities"] == capabilities
        assert (
            dict(
                first.execute(
                    "SELECT table_name, table_role FROM _bioextract.table_info"
                ).fetchall()
            )
            == tables
        )
        assert {
            row[0]
            for row in first.execute(
                "SELECT logical_name FROM _bioextract.source_file"
            ).fetchall()
        } == source_roles
        assert result.tables == tuple(tables)
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_reopened_grouped_selection_preserves_unique_fan_out_and_unmatched(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)

    selection = InterProDatabase.from_duckdb(path).select_groups(
        {"case": ["sp|P12345|TEST", "missing"], "repeat": ["P12345"]}
    )

    assert selection.extract_mapping().select("GroupId", "InputId").to_dicts() == [
        {"GroupId": "case", "InputId": "P12345"},
        {"GroupId": "repeat", "InputId": "P12345"},
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"GroupId": "case", "InputId": "missing"}
    ]


def test_reopened_selection_queries_only_unique_input_ids(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)

    selection = reopened.select_ids(["P12345"])

    def reject_full_mapping(_database: InterProDatabase) -> None:
        raise AssertionError("selection must not extract the full mapping")

    monkeypatch.setattr(InterProDatabase, "extract_mapping", reject_full_mapping)
    assert selection.extract_mapping().height == 1


def test_reopened_empty_selections_preserve_source_typed_outputs(
    tmp_path: Path,
) -> None:
    source = _source(tmp_path)
    path = tmp_path / "interpro.duckdb"
    source.write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)

    for ids in ([], ["", "  "]):
        expected = source.select_ids(ids)
        actual = reopened.select_ids(ids)
        assert actual.extract_mapping().equals(expected.extract_mapping())
        assert actual.extract_unmatched_ids().equals(expected.extract_unmatched_ids())

    expected_grouped = source.select_groups({"empty": ["", "  "]})
    actual_grouped = reopened.select_groups({"empty": ["", "  "]})
    assert actual_grouped.extract_mapping().equals(expected_grouped.extract_mapping())
    assert actual_grouped.extract_unmatched_ids().equals(
        expected_grouped.extract_unmatched_ids()
    )


def test_reopened_full_mapping_preserves_source_deduplication_and_order(
    tmp_path: Path,
) -> None:
    duplicate = "P12345\tIPR000001\tKringle\tPF00051\t10\t80"
    source = _source(tmp_path, mapping_rows=(duplicate, duplicate))
    expected = source.extract_mapping()
    path = tmp_path / "interpro.duckdb"
    source.write_duckdb(path)

    assert InterProDatabase.from_duckdb(path).extract_mapping().equals(expected)


def test_source_handle_connect_reports_missing_capability(tmp_path: Path) -> None:
    with pytest.raises(CapabilityError, match="from_duckdb"):
        _source(tmp_path).connect()


def test_reopened_handle_reports_missing_xml_frame_capability(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)

    with pytest.raises(CapabilityError, match="XML source handle"):
        InterProDatabase.from_duckdb(path).xml_frame("entry")


def test_reopened_handle_reports_missing_publication_source_capability(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)

    with pytest.raises(CapabilityError, match="mapping source"):
        reopened.build_tidy()
    with pytest.raises(CapabilityError, match="mapping source"):
        reopened.write_duckdb(tmp_path / "copy.duckdb")


@pytest.mark.parametrize(
    ("statement", "message"),
    [
        (
            "UPDATE _bioextract.metadata SET value='2' "
            "WHERE key='bioextract.metadata_schema_version'",
            "metadata schema version",
        ),
        (
            "UPDATE _bioextract.metadata SET value='mapping,forged' "
            "WHERE key='bioextract.capabilities'",
            "capability inventory",
        ),
        (
            "UPDATE _bioextract.metadata SET value='forged' "
            "WHERE key='bioextract.source_schema_profile'",
            "source schema profile",
        ),
        (
            "UPDATE _bioextract.table_info SET table_role='forged' "
            "WHERE table_name='mapping'",
            "table role inventory",
        ),
        (
            "UPDATE _bioextract.table_info SET row_count=999 "
            "WHERE table_name='mapping'",
            "row-count drift",
        ),
        (
            "DELETE FROM _bioextract.column_mapping "
            "WHERE table_name='mapping' AND source_column='UniProtId'",
            "column provenance inventory",
        ),
        (
            "UPDATE _bioextract.source_file SET media_type='forged' "
            "WHERE logical_name='protein_to_interpro'",
            "Embedded source inventory",
        ),
    ],
)
def test_from_duckdb_rejects_forged_publication(
    tmp_path: Path,
    statement: str,
    message: str,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(statement)

    with pytest.raises(IntegrityError, match=message):
        InterProDatabase.from_duckdb(path)


def test_atomic_if_exists_preserves_previous_publication(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    source = _source(tmp_path)
    source.write_duckdb(path)
    before = path.read_bytes()

    with pytest.raises(FileExistsError):
        source.write_duckdb(path)
    assert path.read_bytes() == before


@pytest.mark.parametrize(
    "statement",
    [
        "DELETE FROM _bioextract.metadata WHERE key LIKE 'bioextract.release_version%'",
        "UPDATE _bioextract.metadata SET value='caller' "
        "WHERE key='bioextract.release_version_source'",
    ],
)
def test_xml_profile_requires_official_release_metadata(
    tmp_path: Path,
    statement: str,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(statement)

    with pytest.raises(IntegrityError, match="official release metadata"):
        InterProDatabase.from_duckdb(path)


def test_mapping_only_profile_rejects_forged_release_metadata(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path, with_xml=False).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "INSERT INTO _bioextract.metadata VALUES "
            "('bioextract.release_version', '108.0'), "
            "('bioextract.release_version_source', 'caller')"
        )

    with pytest.raises(IntegrityError, match="mapping-only publication"):
        InterProDatabase.from_duckdb(path)


@pytest.mark.parametrize("column", ["interpro_type", "member_db"])
def test_xml_profile_rejects_incomplete_mapping_enrichment(
    tmp_path: Path,
    column: str,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(f'UPDATE mapping SET "{column}"=NULL')

    with pytest.raises(IntegrityError, match="incomplete mapping enrichment"):
        InterProDatabase.from_duckdb(path)


@pytest.mark.parametrize(
    "statement",
    [
        "UPDATE protein_term SET pfam_id='PF99999'",
        "UPDATE term SET pfam_name=''",
        "INSERT INTO term SELECT * FROM term LIMIT 1",
    ],
)
def test_xml_profile_rejects_inconsistent_compact_pfam_relations(
    tmp_path: Path,
    statement: str,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(statement)
        connection.execute(
            "UPDATE _bioextract.table_info SET row_count=(SELECT count(*) FROM term) "
            "WHERE table_name='term'"
        )

    with pytest.raises(IntegrityError, match="compact Pfam relations"):
        InterProDatabase.from_duckdb(path)


@pytest.mark.parametrize("value", ["null", "{}"])
def test_from_duckdb_rejects_wrong_embedded_source_json_shape(
    tmp_path: Path,
    value: str,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? WHERE key='bioextract.sources'",
            [value],
        )

    with pytest.raises(IntegrityError):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_unrecorded_main_schema_view(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("CREATE VIEW forged AS SELECT * FROM mapping")

    with pytest.raises(IntegrityError, match="table inventory"):
        InterProDatabase.from_duckdb(path)
