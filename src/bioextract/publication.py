"""Inspect self-describing bioextract DuckDB publications."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import duckdb

from bioextract._publication import (
    METADATA_SCHEMA_VERSION,
    validate_duckdb_metadata_v3,
)

__all__ = [
    "PublicationDescriptor",
    "PublicationSource",
    "PublicationTable",
    "inspect_duckdb_publication",
]

_INTERNAL_TABLES = {
    "column_mapping",
    "metadata",
    "source_file",
    "table_info",
    "validation_issue",
}
_TABLE_NAME = re.compile(r"^[a-z][a-z0-9_]*$")


@dataclass(frozen=True, slots=True)
class PublicationSource:
    """Describe one physical source recorded by a publication.

    Examples:
        >>> source = PublicationSource(
        ...     logical_name="mapping",
        ...     display_path="raw/mapping.tsv.gz",
        ...     bytes=1024,
        ...     media_type="application/gzip",
        ... )
        >>> source.logical_name
        'mapping'
    """

    logical_name: str
    display_path: str
    bytes: int
    media_type: str
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class PublicationTable:
    """Describe one biological table recorded by a publication.

    Examples:
        >>> table = PublicationTable("compound", "canonical", 201934)
        >>> table.row_count
        201934
    """

    name: str
    role: str
    row_count: int


@dataclass(frozen=True, slots=True)
class PublicationDescriptor:
    """Describe the validated identity and inventory of one DuckDB file.

    The descriptor summarizes discovery-relevant metadata. Complete provenance,
    column mappings, and validation issue details remain authoritative inside
    the inspected DuckDB publication.

    Examples:
        >>> descriptor = PublicationDescriptor(
        ...     path=Path("data.duckdb"),
        ...     metadata_schema_version="3",
        ...     resource_name="example",
        ...     resource_schema_version="example-v1",
        ...     source_schema_profile="example-source-v1",
        ...     package_version="0.1.0",
        ...     generated_at="2026-08-04T00:00:00Z",
        ...     validation_status="passed",
        ...     validation_issue_count=0,
        ... )
        >>> descriptor.resource_name
        'example'
    """

    path: Path
    metadata_schema_version: str
    resource_name: str
    resource_schema_version: str
    source_schema_profile: str
    package_version: str
    generated_at: str
    validation_status: str
    validation_issue_count: int
    source_schema_version: str | None = None
    release_version: str | None = None
    release_version_source: str | None = None
    scope: str | None = None
    sources: tuple[PublicationSource, ...] = ()
    tables: tuple[PublicationTable, ...] = ()


def inspect_duckdb_publication(
    path: os.PathLike[str] | str,
) -> PublicationDescriptor:
    """Validate and describe one bioextract-owned DuckDB publication.

    The file is opened read-only. Inspection validates the current metadata
    schema, embedded/source inventory parity, validation state, biological
    table inventory, and persisted row counts. It does not hash the DuckDB file
    or infer database identity from its path.

    Args:
        path: Existing DuckDB publication to inspect.

    Returns:
        Immutable discovery metadata for the validated publication.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        IsADirectoryError: If ``path`` is not a regular file.
        ValueError: If the file is not a valid current bioextract publication.
        duckdb.Error: If DuckDB cannot open or query the file.

    Examples:
        >>> descriptor = inspect_duckdb_publication(  # doctest: +SKIP
        ...     "tidy/data.duckdb"
        ... )
        >>> descriptor.validation_status  # doctest: +SKIP
        'passed'

    Notes:
        This API owns inspection of one publication only. Resource-tree
        discovery, biofetch manifest parsing, and catalog serialization belong
        to the external publication workflow.
    """
    publication_path = Path(path).resolve(strict=True)
    if not publication_path.is_file():
        raise IsADirectoryError(publication_path)

    connection = duckdb.connect(str(publication_path), read_only=True)
    try:
        connection.execute("PRAGMA database_size").fetchall()
        _validate_internal_table_inventory(connection)

        metadata = dict(
            connection.execute(
                "SELECT key, value FROM _bioextract.metadata ORDER BY key"
            ).fetchall()
        )
        metadata_version = metadata.get("bioextract.metadata_schema_version")
        if metadata_version != METADATA_SCHEMA_VERSION:
            raise ValueError(
                "Unsupported bioextract metadata schema: "
                f"{metadata_version!r}; expected {METADATA_SCHEMA_VERSION!r}"
            )
        validate_duckdb_metadata_v3(connection, metadata)

        sources = tuple(
            PublicationSource(
                logical_name=str(row[0]),
                display_path=str(row[1]),
                bytes=int(row[2]),
                media_type=str(row[3]),
                sha256=None if row[4] is None else str(row[4]),
            )
            for row in connection.execute(
                "SELECT logical_name, display_path, bytes, media_type, sha256 "
                "FROM _bioextract.source_file ORDER BY logical_name"
            ).fetchall()
        )
        tables = _read_and_validate_tables(connection)
        issue_count = int(metadata["bioextract.validation_issue_count"])

        return PublicationDescriptor(
            path=publication_path,
            metadata_schema_version=metadata_version,
            resource_name=metadata["bioextract.resource_name"],
            resource_schema_version=metadata["bioextract.resource_schema_version"],
            source_schema_profile=metadata["bioextract.source_schema_profile"],
            source_schema_version=metadata.get("bioextract.source_schema_version"),
            release_version=metadata.get("bioextract.release_version"),
            release_version_source=metadata.get("bioextract.release_version_source"),
            package_version=metadata["bioextract.package_version"],
            generated_at=metadata["bioextract.generated_at"],
            validation_status=metadata["bioextract.validation_status"],
            validation_issue_count=issue_count,
            scope=metadata.get("bioextract.scope"),
            sources=sources,
            tables=tables,
        )
    finally:
        connection.close()


def _validate_internal_table_inventory(
    connection: duckdb.DuckDBPyConnection,
) -> None:
    rows = connection.execute(
        "SELECT table_name, table_type FROM information_schema.tables "
        "WHERE table_schema='_bioextract' ORDER BY table_name"
    ).fetchall()
    observed = {str(row[0]) for row in rows if row[1] == "BASE TABLE"}
    missing = sorted(_INTERNAL_TABLES - observed)
    unexpected = sorted(observed - _INTERNAL_TABLES)
    if missing or unexpected:
        raise ValueError(
            "Invalid _bioextract table inventory: "
            f"missing={missing}, unexpected={unexpected}"
        )


def _read_and_validate_tables(
    connection: duckdb.DuckDBPyConnection,
) -> tuple[PublicationTable, ...]:
    table_rows = connection.execute(
        "SELECT table_name, table_role, row_count "
        "FROM _bioextract.table_info ORDER BY table_name"
    ).fetchall()
    tables = tuple(
        PublicationTable(str(row[0]), str(row[1]), int(row[2])) for row in table_rows
    )
    declared_names = {table.name for table in tables}
    invalid_names = sorted(
        table.name for table in tables if _TABLE_NAME.fullmatch(table.name) is None
    )
    if invalid_names:
        raise ValueError(f"Publication table names are not query-safe: {invalid_names}")

    observed_names = {
        str(row[0])
        for row in connection.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main' AND table_type='BASE TABLE'"
        ).fetchall()
    }
    if declared_names != observed_names:
        raise ValueError(
            "DuckDB biological table inventory does not match table_info: "
            f"declared={sorted(declared_names)}, actual={sorted(observed_names)}"
        )

    for table in tables:
        count_row = connection.execute(
            f'SELECT count(*) FROM "{table.name}"'
        ).fetchone()
        if count_row is None or int(count_row[0]) != table.row_count:
            raise ValueError(
                f"DuckDB row count does not match table_info for table: {table.name}"
            )
    return tables
