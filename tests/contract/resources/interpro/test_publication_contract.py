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
    xml_entries: str | None = None,
) -> InterProDatabase:
    mapping = tmp_path / "protein2ipr.dat.gz"
    with gzip.open(mapping, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(mapping_rows) + "\n")
    if not with_xml:
        return InterProDatabase.from_mapping_files(protein_to_interpro=mapping)
    xml = tmp_path / "interpro.xml.gz"
    with gzip.open(xml, "wt", encoding="utf-8") as handle:
        handle.write(
            xml_entries
            or """<interprodb>
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

    assert selection.extract_mapping().select("group_id", "input_id").to_dicts() == [
        {"group_id": "case", "input_id": "P12345"},
        {"group_id": "repeat", "input_id": "P12345"},
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"group_id": "case", "input_id": "missing"}
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
            "WHERE table_name='mapping' AND source_column='uniprot_id'",
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


def test_xml_profile_requires_persisted_content_validation_result(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='mapping-v1' "
            "WHERE key='bioextract.interpro_content_validation'"
        )

    with pytest.raises(IntegrityError, match="content-validation result"):
        InterProDatabase.from_duckdb(path)


def test_xml_profile_reopens_multiple_distinct_xrefs_for_one_pfam(
    tmp_path: Path,
) -> None:
    xml = """<interprodb>
<release><dbinfo dbname="INTERPRO" version="108.0"/></release>
<interpro id="IPR000001" type="Domain"><name>First</name><member_list>
<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>
</member_list></interpro>
<interpro id="IPR000002" type="Family"><name>Second</name><member_list>
<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>
</member_list></interpro></interprodb>"""
    source = _source(
        tmp_path,
        mapping_rows=(
            "P12345\tIPR000001\tRaw first name\tPF00051\t10\t80",
            "P12345\tIPR000002\tSecond\tPF00051\t10\t80",
        ),
        xml_entries=xml,
    )
    path = tmp_path / "interpro.duckdb"
    source.write_duckdb(path)

    reopened = InterProDatabase.from_duckdb(path)
    with reopened.connect() as connection:
        assert connection.execute("SELECT count(*) FROM term_xref").fetchone() == (2,)


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


def test_from_duckdb_rejects_unrecorded_provenance_schema_view(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("CREATE VIEW _bioextract.forged AS SELECT * FROM mapping")

    with pytest.raises(IntegrityError, match="provenance table inventory"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_relations_in_additional_user_schema(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE SCHEMA forged; CREATE VIEW forged.mapping AS SELECT * FROM mapping"
        )

    with pytest.raises(IntegrityError, match="unsupported schema"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_normalizes_null_capability_metadata(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE metadata_copy AS SELECT * FROM _bioextract.metadata; "
            "DROP TABLE _bioextract.metadata; "
            "CREATE TABLE _bioextract.metadata AS SELECT * FROM metadata_copy; "
            "UPDATE _bioextract.metadata SET value=NULL "
            "WHERE key='bioextract.capabilities'; "
            "DROP TABLE metadata_copy"
        )

    with pytest.raises(IntegrityError, match="provenance table schema"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_duplicate_metadata_keys(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE metadata_copy AS SELECT * FROM _bioextract.metadata; "
            "DROP TABLE _bioextract.metadata; "
            "CREATE TABLE _bioextract.metadata AS SELECT * FROM metadata_copy; "
            "INSERT INTO _bioextract.metadata VALUES "
            "('bioextract.resource_name', 'interpro'); "
            "DROP TABLE metadata_copy"
        )

    with pytest.raises(IntegrityError, match="duplicate metadata keys"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_duplicate_table_info_keys(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE table_info_copy AS SELECT * FROM _bioextract.table_info; "
            "DROP TABLE _bioextract.table_info; "
            "CREATE TABLE _bioextract.table_info AS SELECT * FROM table_info_copy; "
            "INSERT INTO _bioextract.table_info "
            "SELECT * FROM table_info_copy WHERE table_name='mapping'; "
            "DROP TABLE table_info_copy"
        )

    with pytest.raises(IntegrityError, match="duplicate table-info keys"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_duplicate_column_provenance_keys(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE mapping_copy AS SELECT * FROM _bioextract.column_mapping; "
            "DROP TABLE _bioextract.column_mapping; "
            "CREATE TABLE _bioextract.column_mapping AS SELECT * FROM mapping_copy; "
            "INSERT INTO _bioextract.column_mapping SELECT * FROM mapping_copy LIMIT 1; "
            "DROP TABLE mapping_copy"
        )

    with pytest.raises(IntegrityError, match="duplicate column-provenance keys"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_duplicate_source_file_keys(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE source_copy AS SELECT * FROM _bioextract.source_file; "
            "DROP TABLE _bioextract.source_file; "
            "CREATE TABLE _bioextract.source_file AS SELECT * FROM source_copy; "
            "INSERT INTO _bioextract.source_file SELECT * FROM source_copy LIMIT 1; "
            "DROP TABLE source_copy"
        )

    with pytest.raises(IntegrityError, match="duplicate source-file keys"):
        InterProDatabase.from_duckdb(path)


def test_from_duckdb_rejects_duplicate_validation_issue_keys(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "CREATE TABLE issue_copy AS SELECT * FROM _bioextract.validation_issue; "
            "DROP TABLE _bioextract.validation_issue; "
            "CREATE TABLE _bioextract.validation_issue AS SELECT * FROM issue_copy; "
            "INSERT INTO _bioextract.validation_issue VALUES "
            "(1, 'warning', 'test', 'source', 'mapping', NULL, NULL, NULL, NULL, "
            "NULL, 'test'), "
            "(1, 'warning', 'test', 'source', 'mapping', NULL, NULL, NULL, NULL, "
            "NULL, 'test'); "
            "UPDATE _bioextract.metadata SET value='2' "
            "WHERE key='bioextract.validation_issue_count'; "
            "UPDATE _bioextract.metadata SET value='passed_with_warnings' "
            "WHERE key='bioextract.validation_status'; "
            "DROP TABLE issue_copy"
        )

    with pytest.raises(IntegrityError, match="duplicate validation-issue keys"):
        InterProDatabase.from_duckdb(path)


def test_reopened_handle_rejects_atomically_replaced_publication(
    tmp_path: Path,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)

    replacement_dir = tmp_path / "replacement"
    replacement_dir.mkdir()
    _source(replacement_dir, with_xml=False).write_duckdb(path, if_exists="replace")

    with pytest.raises(IntegrityError, match="was replaced"):
        reopened.connect()


def test_reopened_relative_path_is_bound_to_resolved_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    monkeypatch.chdir(tmp_path)
    reopened = InterProDatabase.from_duckdb("interpro.duckdb")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    with reopened.connect() as connection:
        assert connection.execute("SELECT count(*) FROM mapping").fetchone() == (1,)


def test_reopened_handle_identity_includes_ctime(tmp_path: Path) -> None:
    path = tmp_path / "interpro.duckdb"
    _source(tmp_path).write_duckdb(path)
    reopened = InterProDatabase.from_duckdb(path)
    path.chmod(path.stat().st_mode | 0o100)

    with pytest.raises(IntegrityError, match="was replaced"):
        reopened.connect()


def test_retained_lazy_dataset_rejects_changed_source_before_commit(
    tmp_path: Path,
) -> None:
    source = _source(tmp_path)
    path = tmp_path / "interpro.duckdb"
    source.write_duckdb(path)
    before = path.read_bytes()
    dataset = source.build_tidy()

    mapping_path = source.snapshot.file_protein2ipr
    assert mapping_path is not None
    with gzip.open(mapping_path, "at", encoding="utf-8") as handle:
        handle.write("P12345\tIPR000001\tKringle\tPF00051\t10\t80\n")

    with pytest.raises(IntegrityError, match="source changed"):
        dataset.write_duckdb(path, if_exists="replace")
    assert path.read_bytes() == before


def test_xml_cache_rejects_source_identity_change(tmp_path: Path) -> None:
    source = _source(tmp_path)
    assert source.xml_frame("entry").height == 1
    xml_path = source.snapshot.file_interpro_xml
    assert xml_path is not None
    with gzip.open(xml_path, "at", encoding="utf-8") as handle:
        handle.write("\n")

    with pytest.raises(IntegrityError, match="XML source changed"):
        source.xml_frame("entry")
