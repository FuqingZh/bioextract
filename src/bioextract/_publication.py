from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import duckdb
import polars as pl

METADATA_SCHEMA_VERSION = "3"
_SQL_IDENTIFIER = re.compile(r"^[a-z][a-z0-9_]*$")


@dataclass(frozen=True, slots=True)
class SourceFileRecord:
    """Describe one local source included in publication provenance."""

    logical_name: str
    path: Path
    media_type: str
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class ParquetWriteResult:
    """Describe one successfully committed Parquet publication."""

    path: Path
    resource_name: str
    resource_schema_version: str


@dataclass(frozen=True, slots=True)
class DuckDBWriteResult:
    """Describe one successfully committed multi-relation DuckDB publication."""

    path: Path
    resource_name: str
    resource_schema_version: str
    tables: tuple[str, ...]
    row_counts: Mapping[str, int]
    validation_issue_count: int = 0


@dataclass(frozen=True, slots=True)
class RelationSpec:
    """Bind one lazy relation to its DuckDB table name and semantic role."""

    table_name: str
    frame: pl.LazyFrame
    role: str = "canonical"
    preserve_source_headers: bool = False


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    """Describe one non-fatal source-integrity problem in a publication."""

    severity: str
    issue_code: str
    source_name: str
    relation_name: str
    identifier_namespace: str | None = None
    identifier_value: str | None = None
    referenced_relation: str | None = None
    referenced_identifier: str | None = None
    source_record_number: int | None = None
    message: str = ""


def validate_duckdb_metadata_v3(
    connection: duckdb.DuckDBPyConnection,
    metadata: Mapping[str, str],
) -> None:
    """Validate required v3 keys and embedded/source-table inventory parity."""
    required = {
        "bioextract.resource_name",
        "bioextract.resource_schema_version",
        "bioextract.source_schema_profile",
        "bioextract.package_version",
        "bioextract.generated_at",
        "bioextract.validation_status",
        "bioextract.validation_issue_count",
        "bioextract.sources",
    }
    missing = sorted(required - set(metadata))
    if missing:
        raise ValueError(f"Metadata v3 is missing required keys: {missing}")
    source_rows = connection.execute(
        "SELECT logical_name, display_path, bytes, media_type, sha256 "
        "FROM _bioextract.source_file ORDER BY logical_name"
    ).fetchall()
    embedded_sources = json.loads(metadata["bioextract.sources"])
    table_sources = [
        {
            "logical_name": row[0],
            "path": row[1],
            "bytes": int(row[2]),
            "media_type": row[3],
            **({"sha256": row[4]} if row[4] is not None else {}),
        }
        for row in source_rows
    ]
    if sorted(embedded_sources, key=lambda item: item["logical_name"]) != table_sources:
        raise ValueError("Embedded source inventory does not match source_file")


def write_parquet_publication(
    frame: pl.LazyFrame,
    path: os.PathLike[str] | str,
    *,
    resource_name: str,
    resource_schema_version: str,
    source_schema_profile: str,
    source_schema_version: str | None = None,
    sources: Sequence[SourceFileRecord],
    scope: str | None = None,
    release_version: str | None = None,
    release_version_source: str | None = None,
    if_exists: str = "fail",
    normalize_columns: bool = True,
) -> ParquetWriteResult:
    """Atomically stream one relation to Parquet with embedded provenance."""
    destination = _prepare_destination(path, if_exists=if_exists)
    column_mappings: tuple[tuple[str, str, str, str], ...] = ()
    if normalize_columns:
        frame, column_mappings = _normalize_lazy_columns(frame, table_name="data")
    publication_metadata = _publication_metadata(
        resource_name=resource_name,
        resource_schema_version=resource_schema_version,
        source_schema_profile=source_schema_profile,
        source_schema_version=source_schema_version,
        sources=sources,
        scope=scope,
        release_version=release_version,
        release_version_source=release_version_source,
    )
    publication_metadata["bioextract.column_mapping"] = json.dumps(
        [
            {
                "source_column": source_column,
                "output_column": output_column,
                "reason": reason,
            }
            for _, source_column, output_column, reason in column_mappings
        ],
        separators=(",", ":"),
        sort_keys=True,
    )
    stage = _create_stage_path(destination, suffix=".parquet")
    try:
        frame.sink_parquet(stage, metadata=publication_metadata)
        _validate_parquet_publication(stage, publication_metadata)
        os.replace(stage, destination)
    finally:
        stage.unlink(missing_ok=True)
    return ParquetWriteResult(
        path=destination,
        resource_name=resource_name,
        resource_schema_version=resource_schema_version,
    )


def write_duckdb_publication(
    relations: Sequence[RelationSpec],
    path: os.PathLike[str] | str,
    *,
    resource_name: str,
    resource_schema_version: str,
    source_schema_profile: str,
    source_schema_version: str | None = None,
    sources: Sequence[SourceFileRecord],
    scope: str | None = None,
    release_version: str | None = None,
    release_version_source: str | None = None,
    if_exists: str = "fail",
    column_mappings: Sequence[tuple[str, str, str, str]] = (),
    validation_issues: Sequence[ValidationIssue] = (),
    extra_metadata: Mapping[str, str] | None = None,
) -> DuckDBWriteResult:
    """Atomically publish related lazy frames as one provenance-aware DuckDB."""
    if not relations:
        raise ValueError("At least one relation is required for DuckDB publication")
    for relation in relations:
        _validate_identifier(relation.table_name)

    destination = _prepare_destination(path, if_exists=if_exists)
    stage = _create_stage_path(destination, suffix=".duckdb")
    stage.unlink()
    row_counts: dict[str, int] = {}
    try:
        with tempfile.TemporaryDirectory(prefix="bioextract-relations-") as dir_tmp:
            connection = duckdb.connect(str(stage))
            try:
                _create_metadata_schema(connection)
                mappings = list(column_mappings)
                for index, relation in enumerate(relations):
                    file_relation = Path(dir_tmp) / f"{index:04d}.parquet"
                    frame = relation.frame
                    if relation.preserve_source_headers:
                        frame, relation_mappings = _disambiguate_source_columns(
                            frame,
                            table_name=relation.table_name,
                        )
                        mappings.extend(relation_mappings)
                    else:
                        frame, relation_mappings = _normalize_lazy_columns(
                            frame,
                            table_name=relation.table_name,
                        )
                        mappings.extend(relation_mappings)
                    frame.sink_parquet(file_relation)
                    connection.execute(
                        f'CREATE TABLE "{relation.table_name}" AS '
                        "SELECT * FROM read_parquet(?)",
                        [str(file_relation)],
                    )
                    count_row = connection.execute(
                        f'SELECT count(*) FROM "{relation.table_name}"'
                    ).fetchone()
                    if count_row is None:
                        raise RuntimeError(
                            f"Cannot count published table: {relation.table_name}"
                        )
                    row_counts[relation.table_name] = int(count_row[0])
                publication_metadata = _publication_metadata(
                    resource_name=resource_name,
                    resource_schema_version=resource_schema_version,
                    source_schema_profile=source_schema_profile,
                    source_schema_version=source_schema_version,
                    sources=sources,
                    scope=scope,
                    release_version=release_version,
                    release_version_source=release_version_source,
                    validation_issue_count=len(validation_issues),
                )
                extra = {} if extra_metadata is None else dict(extra_metadata)
                reserved = sorted(set(publication_metadata) & set(extra))
                if reserved:
                    raise ValueError(
                        "extra_metadata cannot replace canonical publication "
                        f"metadata keys: {reserved}"
                    )
                publication_metadata.update(extra)
                _write_duckdb_metadata(
                    connection,
                    metadata=publication_metadata,
                    sources=sources,
                    relations=relations,
                    row_counts=row_counts,
                    column_mappings=mappings,
                    validation_issues=validation_issues,
                )
                connection.execute("CHECKPOINT")
            finally:
                connection.close()
            _validate_duckdb_publication(
                stage,
                relations=relations,
                row_counts=row_counts,
            )
        os.replace(stage, destination)
    finally:
        stage.unlink(missing_ok=True)

    return DuckDBWriteResult(
        path=destination,
        resource_name=resource_name,
        resource_schema_version=resource_schema_version,
        tables=tuple(relation.table_name for relation in relations),
        row_counts=row_counts,
        validation_issue_count=len(validation_issues),
    )


def _publication_metadata(
    *,
    resource_name: str,
    resource_schema_version: str,
    source_schema_profile: str,
    source_schema_version: str | None,
    sources: Sequence[SourceFileRecord],
    scope: str | None,
    release_version: str | None,
    release_version_source: str | None,
    validation_issue_count: int = 0,
) -> dict[str, str]:
    generated_at = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    if not resource_schema_version.strip():
        raise ValueError("resource_schema_version must be non-empty")
    if not source_schema_profile.strip():
        raise ValueError("source_schema_profile must be non-empty")
    metadata = {
        "bioextract.metadata_schema_version": METADATA_SCHEMA_VERSION,
        "bioextract.resource_name": resource_name,
        "bioextract.resource_schema_version": resource_schema_version,
        "bioextract.source_schema_profile": source_schema_profile,
        "bioextract.package_version": _package_version(),
        "bioextract.generated_at": generated_at,
        "bioextract.validation_status": (
            "passed_with_warnings" if validation_issue_count else "passed"
        ),
        "bioextract.validation_issue_count": str(validation_issue_count),
        "bioextract.sources": json.dumps(
            [_source_payload(source) for source in sources],
            separators=(",", ":"),
            sort_keys=True,
        ),
    }
    if scope is not None:
        metadata["bioextract.scope"] = scope
    if source_schema_version is not None:
        if not source_schema_version.strip():
            raise ValueError("source_schema_version must be non-empty when provided")
        metadata["bioextract.source_schema_version"] = source_schema_version
    if release_version is not None:
        if not release_version.strip():
            raise ValueError("release_version must be non-empty when provided")
        if release_version_source is not None and release_version_source not in {
            "caller",
            "official_metadata",
        }:
            raise ValueError(
                "release_version_source must be caller or official_metadata"
            )
        metadata["bioextract.release_version"] = release_version
        metadata["bioextract.release_version_source"] = (
            release_version_source or "caller"
        )
    elif release_version_source is not None:
        raise ValueError("release_version_source requires release_version")
    return metadata


def _create_metadata_schema(connection: duckdb.DuckDBPyConnection) -> None:
    connection.execute("CREATE SCHEMA _bioextract")
    connection.execute(
        """
        CREATE TABLE _bioextract.metadata (
            key VARCHAR PRIMARY KEY,
            value VARCHAR NOT NULL
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE _bioextract.source_file (
            logical_name VARCHAR PRIMARY KEY,
            display_path VARCHAR NOT NULL,
            bytes UBIGINT NOT NULL,
            media_type VARCHAR NOT NULL,
            sha256 VARCHAR
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE _bioextract.table_info (
            table_name VARCHAR PRIMARY KEY,
            table_role VARCHAR NOT NULL,
            row_count UBIGINT NOT NULL
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE _bioextract.column_mapping (
            table_name VARCHAR NOT NULL,
            source_column VARCHAR NOT NULL,
            output_column VARCHAR NOT NULL,
            reason VARCHAR NOT NULL,
            PRIMARY KEY (table_name, source_column)
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE _bioextract.validation_issue (
            issue_id UBIGINT PRIMARY KEY,
            severity VARCHAR NOT NULL,
            issue_code VARCHAR NOT NULL,
            source_name VARCHAR NOT NULL,
            relation_name VARCHAR NOT NULL,
            identifier_namespace VARCHAR,
            identifier_value VARCHAR,
            referenced_relation VARCHAR,
            referenced_identifier VARCHAR,
            source_record_number UBIGINT,
            message VARCHAR NOT NULL
        )
        """
    )


def _write_duckdb_metadata(
    connection: duckdb.DuckDBPyConnection,
    *,
    metadata: Mapping[str, str],
    sources: Sequence[SourceFileRecord],
    relations: Sequence[RelationSpec],
    row_counts: Mapping[str, int],
    column_mappings: Sequence[tuple[str, str, str, str]],
    validation_issues: Sequence[ValidationIssue],
) -> None:
    connection.executemany(
        "INSERT INTO _bioextract.metadata VALUES (?, ?)",
        sorted(metadata.items()),
    )
    connection.executemany(
        "INSERT INTO _bioextract.source_file VALUES (?, ?, ?, ?, ?)",
        [
            (
                source.logical_name,
                str(source.path),
                source.path.stat().st_size,
                source.media_type,
                source.sha256,
            )
            for source in sources
        ],
    )
    connection.executemany(
        "INSERT INTO _bioextract.table_info VALUES (?, ?, ?)",
        [
            (
                relation.table_name,
                relation.role,
                row_counts[relation.table_name],
            )
            for relation in relations
        ],
    )
    if column_mappings:
        connection.executemany(
            "INSERT INTO _bioextract.column_mapping VALUES (?, ?, ?, ?)",
            list(column_mappings),
        )
    if validation_issues:
        connection.executemany(
            "INSERT INTO _bioextract.validation_issue VALUES "
            "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                (
                    issue_id,
                    issue.severity,
                    issue.issue_code,
                    issue.source_name,
                    issue.relation_name,
                    issue.identifier_namespace,
                    issue.identifier_value,
                    issue.referenced_relation,
                    issue.referenced_identifier,
                    issue.source_record_number,
                    issue.message,
                )
                for issue_id, issue in enumerate(validation_issues, start=1)
            ],
        )


def _source_payload(source: SourceFileRecord) -> dict[str, str | int]:
    return {
        "logical_name": source.logical_name,
        "path": str(source.path),
        "bytes": source.path.stat().st_size,
        "media_type": source.media_type,
        **({"sha256": source.sha256} if source.sha256 is not None else {}),
    }


def _validate_parquet_publication(
    path: Path,
    metadata: Mapping[str, str],
) -> None:
    pl.scan_parquet(path).collect_schema()
    connection = duckdb.connect()
    try:
        rows = connection.execute(
            "SELECT CAST(key AS VARCHAR), CAST(value AS VARCHAR) "
            "FROM parquet_kv_metadata(?) "
            "WHERE CAST(key AS VARCHAR) LIKE 'bioextract.%'",
            [str(path)],
        ).fetchall()
    finally:
        connection.close()
    observed = dict(rows)
    missing = sorted(set(metadata) - set(observed))
    if missing:
        raise RuntimeError(f"Parquet publication is missing provenance keys: {missing}")


def _validate_duckdb_publication(
    path: Path,
    *,
    relations: Sequence[RelationSpec],
    row_counts: Mapping[str, int],
) -> None:
    connection = duckdb.connect(str(path), read_only=True)
    try:
        connection.execute("PRAGMA database_size").fetchall()
        metadata_rows = dict(
            connection.execute(
                "SELECT table_name, row_count FROM _bioextract.table_info"
            ).fetchall()
        )
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        issue_row = connection.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone()
        if issue_row is None:
            raise RuntimeError("Cannot count DuckDB validation issues")
        issue_count = int(issue_row[0])
        if (
            metadata.get("bioextract.metadata_schema_version")
            != METADATA_SCHEMA_VERSION
        ):
            raise RuntimeError(
                f"DuckDB publication metadata schema is not v{METADATA_SCHEMA_VERSION}"
            )
        required_metadata = {
            "bioextract.resource_name",
            "bioextract.resource_schema_version",
            "bioextract.source_schema_profile",
            "bioextract.package_version",
            "bioextract.generated_at",
            "bioextract.validation_status",
            "bioextract.validation_issue_count",
            "bioextract.sources",
        }
        missing_metadata = sorted(required_metadata - set(metadata))
        if missing_metadata:
            raise RuntimeError(
                f"DuckDB publication is missing metadata keys: {missing_metadata}"
            )
        source_rows = connection.execute(
            "SELECT logical_name, display_path, bytes, media_type, sha256 "
            "FROM _bioextract.source_file ORDER BY logical_name"
        ).fetchall()
        embedded_sources = json.loads(metadata["bioextract.sources"])
        table_sources = [
            {
                "logical_name": row[0],
                "path": row[1],
                "bytes": int(row[2]),
                "media_type": row[3],
                **({"sha256": row[4]} if row[4] is not None else {}),
            }
            for row in source_rows
        ]
        if (
            sorted(embedded_sources, key=lambda item: item["logical_name"])
            != table_sources
        ):
            raise RuntimeError(
                "DuckDB embedded source inventory does not match source_file"
            )
        if int(metadata.get("bioextract.validation_issue_count", "-1")) != issue_count:
            raise RuntimeError(
                "DuckDB validation issue count does not match publication metadata"
            )
        expected = {
            relation.table_name: row_counts[relation.table_name]
            for relation in relations
        }
        if metadata_rows != expected:
            raise RuntimeError(
                "DuckDB table inventory does not match published row counts"
            )
        for table_name, expected_count in expected.items():
            count_row = connection.execute(
                f'SELECT count(*) FROM "{table_name}"'
            ).fetchone()
            if count_row is None or int(count_row[0]) != expected_count:
                raise RuntimeError(
                    f"DuckDB row-count validation failed for table: {table_name}"
                )
    finally:
        connection.close()


def _normalize_lazy_columns(
    frame: pl.LazyFrame,
    *,
    table_name: str,
) -> tuple[pl.LazyFrame, tuple[tuple[str, str, str, str], ...]]:
    columns = frame.collect_schema().names()
    output_by_source: dict[str, str] = {}
    used: set[str] = set()
    mappings: list[tuple[str, str, str, str]] = []
    for source_column in columns:
        base = _to_snake_case(source_column)
        output_column = base
        index = 2
        while output_column.casefold() in used:
            output_column = f"{base}_{index}"
            index += 1
        used.add(output_column.casefold())
        output_by_source[source_column] = output_column
        if output_column != source_column:
            reason = (
                "identifier_collision"
                if output_column != base
                else "generated_snake_case"
            )
            mappings.append((table_name, source_column, output_column, reason))
    return frame.rename(output_by_source), tuple(mappings)


def _disambiguate_source_columns(
    frame: pl.LazyFrame,
    *,
    table_name: str,
) -> tuple[pl.LazyFrame, tuple[tuple[str, str, str, str], ...]]:
    """Retain official headers except where DuckDB cannot distinguish them."""
    columns = frame.collect_schema().names()
    output_by_source: dict[str, str] = {}
    used: set[str] = set()
    mappings: list[tuple[str, str, str, str]] = []
    for source_column in columns:
        base = source_column if source_column else "column"
        output_column = base
        index = 2
        while output_column.casefold() in used:
            output_column = f"{base}_{index}"
            index += 1
        used.add(output_column.casefold())
        output_by_source[source_column] = output_column
        if output_column != source_column:
            mappings.append(
                (
                    table_name,
                    source_column,
                    output_column,
                    "empty_header"
                    if not source_column
                    else "case_insensitive_collision",
                )
            )
    return frame.rename(output_by_source), tuple(mappings)


def _to_snake_case(value: str) -> str:
    prepared = (
        value.strip()
        .replace("UniProt", "Uniprot")
        .replace("InterPro", "Interpro")
        .replace("WikiPathways", "Wiki_Pathways")
        .replace("EggNOG", "Eggnog")
        .replace("KEGG", "Kegg")
        .replace("STRING", "String")
    )
    prepared = re.sub(r"[^0-9A-Za-z]+", "_", prepared)
    prepared = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", prepared)
    prepared = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", prepared)
    normalized = prepared.strip("_").lower()
    if not normalized:
        raise ValueError(f"Column name cannot be normalized safely: {value!r}")
    if normalized[0].isdigit():
        normalized = f"column_{normalized}"
    return normalized


def _prepare_destination(
    path: os.PathLike[str] | str,
    *,
    if_exists: str,
) -> Path:
    if if_exists not in {"fail", "replace"}:
        raise ValueError("if_exists must be 'fail' or 'replace'")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and if_exists == "fail":
        raise FileExistsError(destination)
    return destination


def _create_stage_path(destination: Path, *, suffix: str) -> Path:
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.",
        suffix=suffix,
        dir=destination.parent,
        delete=False,
    ) as handle:
        return Path(handle.name)


def _validate_identifier(identifier: str) -> None:
    if _SQL_IDENTIFIER.fullmatch(identifier) is None:
        raise ValueError(
            "DuckDB table names must be lowercase snake_case identifiers: "
            f"{identifier!r}"
        )


def _package_version() -> str:
    try:
        return version("bioextract")
    except PackageNotFoundError:
        return "unknown"
