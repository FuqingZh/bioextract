from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import duckdb
import polars as pl
import pytest

import bioextract.publication as publication
from bioextract._publication import (
    RelationSpec,
    SourceFileRecord,
    ValidationIssue,
    write_duckdb_publication,
)
from bioextract.errors import IntegrityError


def _publication(
    tmp_path: Path,
    *,
    validation_issues: tuple[ValidationIssue, ...] = (),
) -> Path:
    source = tmp_path / "source.tsv"
    source.write_text("id\nA\n", encoding="utf-8")
    path = tmp_path / "publication.duckdb"
    write_duckdb_publication(
        (
            RelationSpec("z_relation", pl.DataFrame({"id": [1, 2]}).lazy()),
            RelationSpec("a_relation", pl.DataFrame({"id": [3]}).lazy()),
        ),
        path,
        resource_name="fixture",
        resource_schema_version="fixture-v1",
        source_schema_profile="fixture-source-v1",
        source_schema_version="official-v2",
        sources=(SourceFileRecord("source", source, "text/tab-separated-values"),),
        scope="canonical",
        release_version="2026-08",
        release_version_source="official_metadata",
        column_mappings=(
            ("z_relation", "Z ID", "z_id", "generated_snake_case"),
            ("a_relation", "A ID", "a_id", "generated_snake_case"),
        ),
        validation_issues=validation_issues,
    )
    return path


def test_inspection_preserves_provenance_and_deterministic_order(
    tmp_path: Path,
) -> None:
    issue = ValidationIssue(
        severity="warning",
        issue_code="fixture",
        source_name="source",
        relation_name="z_relation",
        message="fixture warning",
    )
    path = _publication(tmp_path, validation_issues=(issue,))

    result = publication.inspect_publication(path)

    assert result.path == path.resolve()
    assert result.resource_name == "fixture"
    assert result.resource_schema_version == "fixture-v1"
    assert result.source_schema_profile == "fixture-source-v1"
    assert result.source_schema_version == "official-v2"
    assert result.release_version == "2026-08"
    assert result.release_version_source == "official_metadata"
    assert result.scope == "canonical"
    assert result.validation_status == "passed_with_warnings"
    assert result.validation_issue_count == 1
    assert result.table_counts_verified is False
    assert [record.key for record in result.metadata] == sorted(
        record.key for record in result.metadata
    )
    assert [table.table_name for table in result.tables] == [
        "a_relation",
        "z_relation",
    ]
    assert result.source_files[0].display_path == str(tmp_path / "source.tsv")
    assert [mapping.table_name for mapping in result.column_mappings] == [
        "a_relation",
        "z_relation",
    ]
    assert result.validation_issues[0].message == "fixture warning"


def test_default_inspection_opens_read_only_closes_and_does_not_count_domain_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _publication(tmp_path)
    original_connect = duckdb.connect
    observed: dict[str, Any] = {"queries": [], "closed": False}

    class ConnectionProxy:
        def __init__(self, connection: duckdb.DuckDBPyConnection) -> None:
            self.connection = connection

        def execute(self, query: str, *args: Any, **kwargs: Any) -> ConnectionProxy:
            observed["queries"].append(query)
            self.connection.execute(query, *args, **kwargs)
            return self

        def fetchall(self) -> list[Any]:
            return self.connection.fetchall()

        def fetchone(self) -> Any:
            return self.connection.fetchone()

        def close(self) -> None:
            observed["closed"] = True
            self.connection.close()

    def connect(
        path_value: str,
        *,
        read_only: bool = False,
        config: dict[str, Any] | None = None,
    ) -> ConnectionProxy:
        assert read_only is True
        assert config is not None
        assert int(config["threads"]) > 0
        return ConnectionProxy(
            original_connect(path_value, read_only=read_only, config=config)
        )

    monkeypatch.setattr(publication.duckdb, "connect", connect)
    publication.inspect_publication(path)

    assert observed["closed"] is True
    domain_counts = [
        query
        for query in observed["queries"]
        if "count(*)" in query.lower() and "_bioextract" not in query.lower()
    ]
    assert domain_counts == []


def test_count_verification_and_mismatch(tmp_path: Path) -> None:
    path = _publication(tmp_path)
    verified = publication.inspect_publication(path, verify_table_counts=True)
    assert verified.table_counts_verified is True

    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info SET row_count=99 "
            "WHERE table_name='a_relation'"
        )
    with pytest.raises(IntegrityError, match="row count.*a_relation") as caught:
        publication.inspect_publication(path, verify_table_counts=True)
    assert isinstance(caught.value.__cause__, ValueError)


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        (
            "UPDATE _bioextract.metadata SET value='[]' WHERE key='bioextract.sources'",
            "source inventory",
        ),
        (
            "INSERT INTO _bioextract.table_info VALUES "
            "('missing_relation', 'canonical', 0)",
            "does not match biological tables",
        ),
        (
            "UPDATE _bioextract.metadata SET value='2' "
            "WHERE key='bioextract.validation_issue_count'",
            "validation_issue_count",
        ),
    ],
)
def test_inspection_rejects_corrupt_publications(
    tmp_path: Path, corruption: str, message: str
) -> None:
    path = _publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(corruption)
    with pytest.raises(IntegrityError, match=message) as caught:
        publication.inspect_publication(path)
    assert caught.value.__cause__ is not None


def test_inspection_rejects_metadata_only_publication(tmp_path: Path) -> None:
    path = _publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE a_relation")
        connection.execute("DROP TABLE z_relation")
        connection.execute("DELETE FROM _bioextract.table_info")
    with pytest.raises(IntegrityError, match="at least one biological table"):
        publication.inspect_publication(path)


def test_read_only_file_is_inspectable(tmp_path: Path) -> None:
    path = _publication(tmp_path)
    path.chmod(0o444)
    try:
        assert publication.inspect_publication(path).resource_name == "fixture"
    finally:
        path.chmod(0o644)


def test_connection_closes_when_validation_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.source_file")
    original_connect = duckdb.connect
    closed = False

    class ConnectionProxy:
        def __init__(self, connection: duckdb.DuckDBPyConnection) -> None:
            self.connection = connection

        def execute(self, query: str, *args: Any, **kwargs: Any) -> ConnectionProxy:
            self.connection.execute(query, *args, **kwargs)
            return self

        def fetchall(self) -> list[Any]:
            return self.connection.fetchall()

        def close(self) -> None:
            nonlocal closed
            closed = True
            self.connection.close()

    def connect(
        path_value: str,
        *,
        read_only: bool = False,
        config: dict[str, Any] | None = None,
    ) -> ConnectionProxy:
        return ConnectionProxy(
            original_connect(path_value, read_only=read_only, config=config)
        )

    monkeypatch.setattr(publication.duckdb, "connect", connect)
    with pytest.raises(IntegrityError):
        publication.inspect_publication(path)
    assert closed is True


def test_inspection_does_not_need_recorded_source_path(tmp_path: Path) -> None:
    path = _publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        old_path = str(tmp_path / "source.tsv")
        display_path = "../provenance-only/missing-source.tsv"
        connection.execute(
            "UPDATE _bioextract.source_file SET display_path=?", [display_path]
        )
        connection.execute(
            "UPDATE _bioextract.metadata SET value=replace(value, ?, ?) "
            "WHERE key='bioextract.sources'",
            [old_path, display_path],
        )
    result = publication.inspect_publication(path)
    assert result.source_files[0].display_path == display_path
    assert not os.path.exists(tmp_path / "provenance-only" / "missing-source.tsv")
