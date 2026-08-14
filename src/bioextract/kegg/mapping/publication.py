from __future__ import annotations

from pathlib import Path
from typing import Literal

import duckdb
import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    validate_duckdb_metadata_v2,
)
from bioextract.errors import IntegrityError

from .constant import (
    CAPABILITY_NAMES,
    SCHEMA_VERSION,
    SOURCE_SCHEMA_PROFILE,
    TABLE_SCHEMAS,
)
from .source import MappingSnapshot

_TABLE_ROLES = {
    "organism": "canonical",
    "gene_annotation": "canonical",
    "ko_annotation": "canonical",
}


def write_mapping_publication(
    snapshot: MappingSnapshot,
    path: Path,
    *,
    if_exists: Literal["fail", "replace"],
) -> DuckDBWriteResult:
    """Publish one KEGG mapping build with native DuckDB scans and SQL."""
    from ._native import write_native_mapping_publication

    return write_native_mapping_publication(snapshot, path, if_exists=if_exists)


def open_mapping_publication(path: Path) -> MappingSnapshot:
    metadata, capabilities, members = validate_mapping_publication(path)
    return MappingSnapshot(
        mode="publication",
        publication_path=path,
        publication_capabilities=capabilities,
        publication_members=members,
        release_version=metadata.get("bioextract.release_version"),
    )


def validate_mapping_publication(
    path: Path,
) -> tuple[dict[str, str], dict[str, bool], tuple[str, ...]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        connection = duckdb.connect(str(path), read_only=True)
    except duckdb.Error as error:
        raise IntegrityError(f"Cannot open KEGG mapping publication: {path}") from error
    try:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        try:
            validate_duckdb_metadata_v2(connection, metadata)
        except (duckdb.Error, ValueError) as error:
            raise IntegrityError(str(error)) from error
        expected_identity = {
            "bioextract.resource_name": "kegg",
            "bioextract.resource_schema_version": SCHEMA_VERSION,
            "bioextract.source_schema_profile": SOURCE_SCHEMA_PROFILE,
            "bioextract.scope": "mapping",
        }
        for key, expected in expected_identity.items():
            if metadata.get(key) != expected:
                raise IntegrityError(
                    f"KEGG mapping publication identity mismatch for {key}"
                )
        capability_keys = {
            key: value
            for key, value in metadata.items()
            if key.startswith("bioextract.capability.")
        }
        expected_keys = {f"bioextract.capability.{name}" for name in CAPABILITY_NAMES}
        if set(capability_keys) != expected_keys:
            raise IntegrityError("KEGG mapping capability inventory is unsupported")
        if any(value not in {"true", "false"} for value in capability_keys.values()):
            raise IntegrityError("KEGG mapping capabilities must be true or false")
        capabilities = {
            key.removeprefix("bioextract.capability."): value == "true"
            for key, value in capability_keys.items()
        }
        physical = {
            str(row[0])
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema='main' AND table_type='BASE TABLE'"
            ).fetchall()
        }
        table_info = {
            str(row[0]): (str(row[1]), int(row[2]))
            for row in connection.execute(
                "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
            ).fetchall()
        }
        if physical != set(TABLE_SCHEMAS) or set(table_info) != set(TABLE_SCHEMAS):
            raise IntegrityError("KEGG mapping table inventory is unsupported")
        for table_name, schema in TABLE_SCHEMAS.items():
            observed = tuple(
                (str(row[0]), str(row[1]))
                for row in connection.execute(f'DESCRIBE "{table_name}"').fetchall()
            )
            expected = tuple(
                (name, _duckdb_type(dtype)) for name, dtype in schema.items()
            )
            if observed != expected:
                raise IntegrityError(
                    f"KEGG mapping table schema is unsupported: {table_name}"
                )
            count = _scalar_int(connection, f'SELECT count(*) FROM "{table_name}"')
            if table_info[table_name] != (_TABLE_ROLES[table_name], count):
                raise IntegrityError(f"KEGG mapping row-count drift: {table_name}")
        members = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT organism_code FROM organism ORDER BY organism_code"
            ).fetchall()
        )
        if not members or len(members) != len(set(members)):
            raise IntegrityError(
                "KEGG mapping organism members are empty or duplicated"
            )
        mismatch = _scalar_int(
            connection,
            "SELECT count(*) FROM gene_annotation "
            "WHERE organism_code IS NULL OR kegg_gene_id IS NULL "
            "OR NOT starts_with(kegg_gene_id, organism_code || ':')",
        )
        if mismatch:
            raise IntegrityError("KEGG gene rows cross organism boundaries")
        _validate_capability_nulls(connection, capabilities)
    except duckdb.Error as error:
        raise IntegrityError(f"Invalid KEGG mapping publication: {path}") from error
    finally:
        connection.close()
    return metadata, capabilities, members


def _validate_capability_nulls(
    connection: duckdb.DuckDBPyConnection, capabilities: dict[str, bool]
) -> None:
    gene_columns = {
        "gene_list": "gene_aliases",
        "uniprot_conversion": "uniprot_mappings",
        "ncbi_gene_conversion": "ncbi_gene_mappings",
        "gene_ko": "ko_mappings",
        "gene_pathway": "pathway_mappings",
    }
    for capability, column in gene_columns.items():
        operator = "IS NULL" if capabilities[capability] else "IS NOT NULL"
        count = _scalar_int(
            connection,
            f'SELECT count(*) FROM gene_annotation WHERE "{column}" {operator}',
        )
        if count:
            raise IntegrityError(
                f"KEGG mapping capability/null invariant failed: {capability}"
            )
    ko_operator = "IS NULL" if capabilities["ko_pathway"] else "IS NOT NULL"
    if _scalar_int(
        connection,
        f"SELECT count(*) FROM ko_annotation WHERE pathway_mappings {ko_operator}",
    ):
        raise IntegrityError(
            "KEGG mapping capability/null invariant failed: ko_pathway"
        )


def _scalar_int(connection: duckdb.DuckDBPyConnection, query: str) -> int:
    row = connection.execute(query).fetchone()
    if row is None:
        raise IntegrityError(f"KEGG mapping validation query returned no row: {query}")
    return int(row[0])


def _duckdb_type(dtype: object) -> str:
    if dtype == pl.String:
        return "VARCHAR"
    if isinstance(dtype, pl.List):
        return f"{_duckdb_type(dtype.inner)}[]"
    if isinstance(dtype, pl.Struct):
        fields = ", ".join(
            f"{field.name} {_duckdb_type(field.dtype)}" for field in dtype.fields
        )
        return f"STRUCT({fields})"
    raise AssertionError(f"Unsupported KEGG mapping dtype: {dtype}")
