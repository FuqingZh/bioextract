from __future__ import annotations

import json
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract._publication import validate_duckdb_metadata_v3
from bioextract._tidy import TidyAsset, TidyDataset, TidySource


def _dataset(tmp_path: Path, *, relation_count: int = 1) -> TidyDataset:
    source = tmp_path / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    frames = {
        f"relation_{index}": pl.DataFrame({"id": [f"T{index}"]}).lazy()
        for index in range(relation_count)
    }
    return TidyDataset(
        frames=frames,
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=tuple(
            TidyAsset(
                path=f"relation_{index}.parquet",
                kind="canonical",
                frame_name=f"relation_{index}",
            )
            for index in range(relation_count)
        ),
        resource_name="example",
        release_version="2026-07-29",
    )


def test_parquet_publication_embeds_provenance_without_sidecar(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.parquet"
    result = _dataset(tmp_path).write_parquet(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    metadata = dict(
        duckdb.connect()
        .execute(
            "SELECT key, value FROM parquet_kv_metadata(?) "
            "WHERE CAST(key AS VARCHAR) LIKE 'bioextract.%'",
            [str(path)],
        )
        .fetchall()
    )
    assert metadata[b"bioextract.resource_name"] == b"example"
    assert metadata[b"bioextract.resource_schema_version"] == b"example-v1"
    assert metadata[b"bioextract.release_version"] == b"2026-07-29"


@pytest.mark.parametrize(
    ("release_version", "release_version_source", "message"),
    [
        (" ", None, "release_version must be non-empty"),
        ("2026_01", "filename", "caller or official_metadata"),
        (None, "caller", "requires release_version"),
    ],
)
def test_publication_rejects_invalid_release_provenance(
    tmp_path: Path,
    release_version: str | None,
    release_version_source: str | None,
    message: str,
) -> None:
    dataset = _dataset(tmp_path)
    dataset.release_version = release_version
    dataset.release_version_source = release_version_source
    with pytest.raises(ValueError, match=message):
        dataset.write_parquet(tmp_path / "invalid-release.parquet")


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        (
            "DELETE FROM _bioextract.metadata "
            "WHERE key='bioextract.release_version_source'",
            "must occur together",
        ),
        (
            "UPDATE _bioextract.metadata SET value='filename' "
            "WHERE key='bioextract.release_version_source'",
            "caller or official_metadata",
        ),
        (
            "UPDATE _bioextract.metadata SET value=' ' "
            "WHERE key='bioextract.release_version'",
            "release_version must be non-empty",
        ),
    ],
)
def test_v3_reader_rejects_invalid_release_provenance(
    tmp_path: Path, corruption: str, message: str
) -> None:
    path = tmp_path / "release.duckdb"
    _dataset(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(corruption)
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        with pytest.raises(ValueError, match=message):
            validate_duckdb_metadata_v3(connection, metadata)


def test_v3_reader_rejects_release_source_without_release(tmp_path: Path) -> None:
    path = tmp_path / "source-only.duckdb"
    dataset = _dataset(tmp_path)
    dataset.release_version = None
    dataset.write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "INSERT INTO _bioextract.metadata VALUES "
            "('bioextract.release_version_source', 'caller')"
        )
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        with pytest.raises(ValueError, match="must occur together"):
            validate_duckdb_metadata_v3(connection, metadata)


@pytest.mark.parametrize(
    ("status", "count", "insert_issue", "message"),
    [
        ("warning", "0", False, "validation_status"),
        ("passed_with_warnings", "0", False, "does not match"),
        ("passed", "1", False, "does not match validation_issue"),
        ("passed", "1", True, "does not match validation_issue_count"),
    ],
)
def test_v3_reader_validates_status_and_issue_count_parity(
    tmp_path: Path,
    status: str,
    count: str,
    insert_issue: bool,
    message: str,
) -> None:
    path = tmp_path / "validation-state.duckdb"
    _dataset(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? "
            "WHERE key='bioextract.validation_status'",
            [status],
        )
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? "
            "WHERE key='bioextract.validation_issue_count'",
            [count],
        )
        if insert_issue:
            connection.execute(
                "INSERT INTO _bioextract.validation_issue "
                "(issue_id, severity, issue_code, source_name, relation_name, message) "
                "VALUES (1, 'warning', 'fixture', 'source', 'relation_0', 'fixture')"
            )
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        with pytest.raises(ValueError, match=message):
            validate_duckdb_metadata_v3(connection, metadata)


def test_v3_reader_reports_missing_validation_issue_table(tmp_path: Path) -> None:
    path = tmp_path / "missing-validation-table.duckdb"
    _dataset(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        connection.execute("DROP TABLE _bioextract.validation_issue")
        with pytest.raises(ValueError, match="requires _bioextract.validation_issue"):
            validate_duckdb_metadata_v3(connection, metadata)


@pytest.mark.parametrize("container", ["parquet", "duckdb"])
@pytest.mark.parametrize("logical_name", ["", "   "])
def test_tidy_source_logical_name_must_be_nonempty_after_normalization(
    tmp_path: Path, container: str, logical_name: str
) -> None:
    dataset = _dataset(tmp_path)
    source = dataset._sources[0]  # pyright: ignore[reportPrivateUsage]
    dataset.source = TidySource(
        logical_name, source.path, source.media_type, source.sha256
    )

    with pytest.raises(ValueError, match="logical_name must be non-empty"):
        if container == "parquet":
            dataset.write_parquet(tmp_path / "invalid.parquet")
        else:
            dataset.write_duckdb(tmp_path / "invalid.duckdb")


@pytest.mark.parametrize("container", ["parquet", "duckdb"])
def test_tidy_source_logical_name_is_normalized_before_publication(
    tmp_path: Path, container: str
) -> None:
    dataset = _dataset(tmp_path)
    source = dataset._sources[0]  # pyright: ignore[reportPrivateUsage]
    dataset.source = TidySource(
        " source ", source.path, source.media_type, source.sha256
    )
    path = tmp_path / f"normalized.{container}"
    if container == "parquet":
        dataset.write_parquet(path)
        with duckdb.connect() as connection:
            value = connection.execute(
                "SELECT decode(value) FROM parquet_kv_metadata(?) "
                "WHERE CAST(key AS VARCHAR)='bioextract.sources'",
                [str(path)],
            ).fetchone()
        assert value is not None
        assert json.loads(value[0])[0]["logical_name"] == "source"
    else:
        dataset.write_duckdb(path)
        with duckdb.connect(str(path), read_only=True) as connection:
            assert connection.execute(
                "SELECT logical_name FROM _bioextract.source_file"
            ).fetchone() == ("source",)


@pytest.mark.parametrize("container", ["parquet", "duckdb"])
def test_tidy_source_rejects_normalized_logical_name_collisions(
    tmp_path: Path, container: str
) -> None:
    dataset = _dataset(tmp_path)
    source = dataset._sources[0]  # pyright: ignore[reportPrivateUsage]
    dataset.source = (
        TidySource("source", source.path, source.media_type),
        TidySource(" source ", source.path, source.media_type),
    )

    with pytest.raises(ValueError, match="unique after normalization"):
        if container == "parquet":
            dataset.write_parquet(tmp_path / "collision.parquet")
        else:
            dataset.write_duckdb(tmp_path / "collision.duckdb")


def test_duckdb_publication_has_internal_provenance_schema(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.duckdb"
    result = _dataset(tmp_path, relation_count=2).write_duckdb(path)

    assert result.tables == ("relation_0", "relation_1")
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata_tables = {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = '_bioextract'"
            ).fetchall()
        }
        assert metadata_tables == {
            "column_mapping",
            "metadata",
            "source_file",
            "table_info",
            "validation_issue",
        }
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key = 'bioextract.metadata_schema_version'"
        ).fetchone() == ("3",)
        assert connection.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone() == (0,)
        resource_row = connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key = 'bioextract.resource_name'"
        ).fetchone()
        assert resource_row is not None
        assert resource_row[0] == "example"
        assert connection.execute(
            "SELECT table_name, table_role, row_count "
            "FROM _bioextract.table_info ORDER BY table_name"
        ).fetchall() == [
            ("relation_0", "canonical", 1),
            ("relation_1", "canonical", 1),
        ]
        mapping_row = connection.execute(
            "SELECT count(*) FROM _bioextract.column_mapping"
        ).fetchone()
        assert mapping_row is not None
        assert mapping_row[0] == 0
