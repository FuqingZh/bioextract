"""Read-only inspection of one explicit bioextract DuckDB publication."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import duckdb

from bioextract._publication import (
    METADATA_SCHEMA_VERSION,
    validate_duckdb_metadata_v1,
)
from bioextract.errors import IntegrityError

__all__ = [
    "PublicationColumnMapping",
    "PublicationInspection",
    "PublicationMetadata",
    "PublicationSourceFile",
    "PublicationTable",
    "PublicationValidationIssue",
    "inspect_publication",
]


@dataclass(frozen=True, slots=True)
class PublicationMetadata:
    """Preserve one embedded publication metadata entry.

    Examples:
        >>> PublicationMetadata("bioextract.metadata_schema_version", "1").value
        '1'
    """

    key: str
    value: str


@dataclass(frozen=True, slots=True)
class PublicationTable:
    """Describe one biological table recorded by the publication.

    Examples:
        >>> PublicationTable("term", "canonical", 12).table_name
        'term'
    """

    table_name: str
    table_role: str
    row_count: int


@dataclass(frozen=True, slots=True)
class PublicationSourceFile:
    """Preserve one source-file provenance record.

    Examples:
        >>> PublicationSourceFile("terms", "inputs/terms.obo", 42, "text/obo", None).bytes
        42
    """

    logical_name: str
    display_path: str
    bytes: int
    media_type: str
    sha256: str | None


@dataclass(frozen=True, slots=True)
class PublicationColumnMapping:
    """Describe one published source-to-output column mapping.

    Examples:
        >>> PublicationColumnMapping("term", "Term ID", "term_id", "generated_snake_case").output_column
        'term_id'
    """

    table_name: str
    source_column: str
    output_column: str
    reason: str


@dataclass(frozen=True, slots=True)
class PublicationValidationIssue:
    """Preserve one non-fatal validation issue recorded by the writer.

    Examples:
        >>> issue = PublicationValidationIssue(1, "warning", "missing", "terms", "term", None, None, None, None, None, "Missing parent")
        >>> issue.issue_code
        'missing'
    """

    issue_id: int
    severity: str
    issue_code: str
    source_name: str
    relation_name: str
    identifier_namespace: str | None
    identifier_value: str | None
    referenced_relation: str | None
    referenced_identifier: str | None
    source_record_number: int | None
    message: str


@dataclass(frozen=True, slots=True)
class PublicationInspection:
    """Validated identity, provenance, and inventory for one publication.

    ``table_counts_verified`` distinguishes structural inspection from an
    explicit scan of every biological table.

    Examples:
        >>> inspection = inspect_publication("publication.duckdb")  # doctest: +SKIP
        >>> inspection.resource_name  # doctest: +SKIP
        'go'
    """

    path: Path
    resource_name: str
    resource_schema_version: str
    source_schema_profile: str
    source_schema_version: str | None
    release_version: str | None
    release_version_source: str | None
    metadata_schema_version: str
    package_version: str
    generated_at: str
    scope: str | None
    validation_status: str
    validation_issue_count: int
    table_counts_verified: bool
    metadata: tuple[PublicationMetadata, ...]
    tables: tuple[PublicationTable, ...]
    source_files: tuple[PublicationSourceFile, ...]
    column_mappings: tuple[PublicationColumnMapping, ...]
    validation_issues: tuple[PublicationValidationIssue, ...]


def inspect_publication(
    path: os.PathLike[str] | str,
    *,
    verify_table_counts: bool = False,
) -> PublicationInspection:
    """Validate and describe one caller-supplied bioextract DuckDB file.

    The file is opened read-only and its connection is closed on both success
    and failure. Default inspection validates publication structure and
    provenance without scanning biological rows. Set ``verify_table_counts``
    to opt into ``count(*)`` validation for every biological table.

    Args:
        path: Exact local DuckDB publication to inspect.
        verify_table_counts: Whether to verify recorded biological row counts.

    Returns:
        Immutable, deterministically ordered publication records.

    Raises:
        TypeError: If ``verify_table_counts`` is not exactly a boolean.
        IntegrityError: If the file cannot be opened or fails validation.

    Examples:
        >>> result = inspect_publication("publication.duckdb")  # doctest: +SKIP
        >>> result.table_counts_verified
        False

    Notes:
        This API inspects exactly one path. It does not discover publications,
        infer identity from filenames, or expose the underlying connection.
    """
    if not isinstance(  # pyright: ignore[reportUnnecessaryIsInstance]
        verify_table_counts, bool
    ):
        raise TypeError("verify_table_counts must be a bool")

    publication_path: Path | None = None
    connection: duckdb.DuckDBPyConnection | None = None
    try:
        publication_path = Path(path).resolve()
        connection = duckdb.connect(str(publication_path), read_only=True)
        metadata_rows = tuple(
            PublicationMetadata(str(key), str(value))
            for key, value in connection.execute(
                "SELECT key, value FROM _bioextract.metadata ORDER BY key"
            ).fetchall()
        )
        metadata_values = {record.key: record.value for record in metadata_rows}
        if (
            metadata_values.get("bioextract.metadata_schema_version")
            != METADATA_SCHEMA_VERSION
        ):
            raise ValueError(
                f"Publication metadata schema must be v{METADATA_SCHEMA_VERSION}"
            )
        validate_duckdb_metadata_v1(connection, metadata_values)
        _validate_identity_values(metadata_values)

        tables = tuple(
            PublicationTable(str(name), str(role), int(count))
            for name, role, count in connection.execute(
                "SELECT table_name, table_role, row_count "
                "FROM _bioextract.table_info ORDER BY table_name"
            ).fetchall()
        )
        _validate_table_inventory(connection, tables)
        if verify_table_counts:
            _verify_table_counts(connection, tables)

        source_files = tuple(
            PublicationSourceFile(
                str(logical_name),
                str(display_path),
                int(byte_count),
                str(media_type),
                None if sha256 is None else str(sha256),
            )
            for logical_name, display_path, byte_count, media_type, sha256 in (
                connection.execute(
                    "SELECT logical_name, display_path, bytes, media_type, sha256 "
                    "FROM _bioextract.source_file ORDER BY logical_name"
                ).fetchall()
            )
        )
        column_mappings = tuple(
            PublicationColumnMapping(
                str(table_name), str(source), str(output), str(reason)
            )
            for table_name, source, output, reason in connection.execute(
                "SELECT table_name, source_column, output_column, reason "
                "FROM _bioextract.column_mapping "
                "ORDER BY table_name, source_column, output_column, reason"
            ).fetchall()
        )
        validation_issues = tuple(
            PublicationValidationIssue(
                int(row[0]),
                str(row[1]),
                str(row[2]),
                str(row[3]),
                str(row[4]),
                None if row[5] is None else str(row[5]),
                None if row[6] is None else str(row[6]),
                None if row[7] is None else str(row[7]),
                None if row[8] is None else str(row[8]),
                None if row[9] is None else int(row[9]),
                str(row[10]),
            )
            for row in connection.execute(
                "SELECT issue_id, severity, issue_code, source_name, relation_name, "
                "identifier_namespace, identifier_value, referenced_relation, "
                "referenced_identifier, source_record_number, message "
                "FROM _bioextract.validation_issue "
                "ORDER BY issue_id, severity, issue_code, source_name, relation_name"
            ).fetchall()
        )
        return PublicationInspection(
            path=publication_path,
            resource_name=metadata_values["bioextract.resource_name"],
            resource_schema_version=metadata_values[
                "bioextract.resource_schema_version"
            ],
            source_schema_profile=metadata_values["bioextract.source_schema_profile"],
            source_schema_version=metadata_values.get(
                "bioextract.source_schema_version"
            ),
            release_version=metadata_values.get("bioextract.release_version"),
            release_version_source=metadata_values.get(
                "bioextract.release_version_source"
            ),
            metadata_schema_version=metadata_values[
                "bioextract.metadata_schema_version"
            ],
            package_version=metadata_values["bioextract.package_version"],
            generated_at=metadata_values["bioextract.generated_at"],
            scope=metadata_values.get("bioextract.scope"),
            validation_status=metadata_values["bioextract.validation_status"],
            validation_issue_count=int(
                metadata_values["bioextract.validation_issue_count"]
            ),
            table_counts_verified=verify_table_counts,
            metadata=metadata_rows,
            tables=tables,
            source_files=source_files,
            column_mappings=column_mappings,
            validation_issues=validation_issues,
        )
    except Exception as error:
        path_context = publication_path if publication_path is not None else repr(path)
        raise IntegrityError(
            f"Cannot inspect bioextract publication at {path_context}: {error}"
        ) from error
    finally:
        if connection is not None:
            connection.close()


def _validate_identity_values(metadata: dict[str, str]) -> None:
    for key in (
        "bioextract.resource_name",
        "bioextract.resource_schema_version",
        "bioextract.source_schema_profile",
        "bioextract.package_version",
        "bioextract.generated_at",
    ):
        if not metadata[key].strip():
            raise ValueError(f"Metadata v1 value must be non-empty: {key}")
    source_schema_version = metadata.get("bioextract.source_schema_version")
    if source_schema_version is not None and not source_schema_version.strip():
        raise ValueError("Metadata v1 source_schema_version must be non-empty")


def _validate_table_inventory(
    connection: duckdb.DuckDBPyConnection,
    tables: tuple[PublicationTable, ...],
) -> None:
    recorded = {table.table_name for table in tables}
    biological = {
        str(row[0])
        for row in connection.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main' AND table_type='BASE TABLE'"
        ).fetchall()
    }
    if not recorded or not biological:
        raise ValueError("Publication requires at least one biological table")
    if recorded != biological:
        raise ValueError(
            "Publication table_info does not match biological tables in main: "
            f"recorded={sorted(recorded)}, actual={sorted(biological)}"
        )


def _verify_table_counts(
    connection: duckdb.DuckDBPyConnection,
    tables: tuple[PublicationTable, ...],
) -> None:
    for table in tables:
        escaped_name = table.table_name.replace('"', '""')
        row = connection.execute(f'SELECT count(*) FROM "{escaped_name}"').fetchone()
        if row is None or int(row[0]) != table.row_count:
            actual = None if row is None else int(row[0])
            raise ValueError(
                "Publication row count does not match table_info: "
                f"table={table.table_name!r}, recorded={table.row_count}, "
                f"actual={actual}"
            )
