from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.publication import (
    PublicationDescriptor,
    PublicationSource,
    PublicationTable,
    inspect_duckdb_publication,
)


def _write_publication(tmp_path: Path) -> Path:
    source = tmp_path / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    path = tmp_path / "data.duckdb"
    TidyDataset(
        frames={"term": pl.DataFrame({"id": ["T1", "T2"]}).lazy()},
        source=TidySource(
            "source",
            source,
            "text/tab-separated-values",
            "a" * 64,
        ),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        source_schema_version="upstream-v2",
        build_id_prefix="example",
        assets=(TidyAsset("term.parquet", "canonical", "term"),),
        resource_name="example",
        release_version="2026-08-04",
        release_version_source="official_metadata",
    ).write_duckdb(path)
    return path


def test_inspect_duckdb_publication_returns_validated_descriptor(
    tmp_path: Path,
) -> None:
    path = _write_publication(tmp_path)
    modified_before = path.stat().st_mtime_ns

    descriptor = inspect_duckdb_publication(path)

    assert descriptor == PublicationDescriptor(
        path=path.resolve(),
        metadata_schema_version="3",
        resource_name="example",
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        source_schema_version="upstream-v2",
        release_version="2026-08-04",
        release_version_source="official_metadata",
        package_version="0.1.0",
        generated_at=descriptor.generated_at,
        validation_status="passed",
        validation_issue_count=0,
        sources=(
            PublicationSource(
                logical_name="source",
                display_path=str(tmp_path / "source.tsv"),
                bytes=6,
                media_type="text/tab-separated-values",
                sha256="a" * 64,
            ),
        ),
        tables=(PublicationTable("term", "canonical", 2),),
    )
    assert path.stat().st_mtime_ns == modified_before


def test_inspection_rejects_non_publication_file(tmp_path: Path) -> None:
    path = tmp_path / "plain.duckdb"
    with duckdb.connect(str(path)) as connection:
        connection.execute("CREATE TABLE data (id INTEGER)")

    with pytest.raises(ValueError, match="Invalid _bioextract table inventory"):
        inspect_duckdb_publication(path)


def test_inspection_rejects_biological_table_inventory_drift(
    tmp_path: Path,
) -> None:
    path = _write_publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("CREATE TABLE undeclared (id INTEGER)")

    with pytest.raises(ValueError, match="does not match table_info"):
        inspect_duckdb_publication(path)


def test_inspection_rejects_persisted_row_count_drift(tmp_path: Path) -> None:
    path = _write_publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info SET row_count=999 WHERE table_name='term'"
        )

    with pytest.raises(ValueError, match="row count does not match"):
        inspect_duckdb_publication(path)


def test_inspection_rejects_unknown_metadata_schema(tmp_path: Path) -> None:
    path = _write_publication(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='999' "
            "WHERE key='bioextract.metadata_schema_version'"
        )

    with pytest.raises(ValueError, match="Unsupported bioextract metadata schema"):
        inspect_duckdb_publication(path)


def test_inspection_requires_regular_existing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        inspect_duckdb_publication(tmp_path / "missing.duckdb")
    with pytest.raises(IsADirectoryError):
        inspect_duckdb_publication(tmp_path)


def test_publication_module_exports_only_stable_inspection_surface() -> None:
    import bioextract.publication as publication

    assert publication.__all__ == [
        "PublicationDescriptor",
        "PublicationSource",
        "PublicationTable",
        "inspect_duckdb_publication",
    ]
