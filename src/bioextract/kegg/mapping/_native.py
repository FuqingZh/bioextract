"""DuckDB-native ingestion for the KEGG organism mapping publication.

The module is deliberately private.  It owns the set-oriented ingestion path
used by the publisher, while source-backed lazy relations continue to use the
replayable Polars adapter in :mod:`bioextract.kegg.mapping.query`.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

import duckdb
import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    _connect_publication,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    _create_metadata_schema,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    _create_stage_path,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    _publication_metadata,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    _validate_duckdb_publication,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    _write_duckdb_metadata,  # pyright: ignore[reportPrivateUsage]  # shared publication boundary
    preflight_publication_destination,
    require_package_version,
)
from bioextract.errors import IntegrityError

from .constant import (
    MEDIA_TYPE_TSV,
    SCHEMA_VERSION,
    SOURCE_SCHEMA_PROFILE,
    TABLE_SCHEMAS,
)
from .source import MappingSnapshot, resolve_organism_work, source_capabilities

_REJECTS_LIMIT = 1024
_AGGREGATION_PARTITIONS = 8
_THREAD_CAP = 8
_TABLE_ROLES = dict.fromkeys(TABLE_SCHEMAS, "canonical")
_ROLE_COLUMNS: dict[str, tuple[str, ...]] = {
    "gene_list": ("gene_id", "gene_type", "genomic_position", "display"),
    "uniprot_conversion": ("xref", "gene_id"),
    "ncbi_gene_conversion": ("xref", "gene_id"),
    "gene_ko": ("gene_id", "ko_raw"),
    "gene_pathway": ("gene_id", "pathway_raw"),
    "organism_list": ("genome_id", "organism_code", "organism_name", "taxonomy"),
    "ko_pathway": ("ko_raw", "pathway_raw"),
}
_ORGANISM_ROLES = (
    "gene_list",
    "uniprot_conversion",
    "ncbi_gene_conversion",
    "gene_ko",
    "gene_pathway",
)
_NORMALIZED_GENE_DEFINITIONS = {
    "gene_list": (
        "_kegg_gene_list",
        "gene_id VARCHAR, gene_type VARCHAR, genomic_position VARCHAR, display VARCHAR",
        "SELECT DISTINCT i.organism_code, i.source_path, trim(r.gene_id) AS gene_id, "
        "nullif(trim(r.gene_type), '') AS gene_type, "
        "nullif(trim(r.genomic_position), '') AS genomic_position, trim(r.display) AS display",
    ),
    "uniprot_conversion": (
        "_kegg_uniprot",
        "xref VARCHAR, gene_id VARCHAR",
        "SELECT DISTINCT i.organism_code, i.source_path, trim(r.xref) AS xref, trim(r.gene_id) AS gene_id",
    ),
    "ncbi_gene_conversion": (
        "_kegg_ncbi",
        "xref VARCHAR, gene_id VARCHAR",
        "SELECT DISTINCT i.organism_code, i.source_path, trim(r.xref) AS xref, trim(r.gene_id) AS gene_id",
    ),
    "gene_ko": (
        "_kegg_gene_ko",
        "gene_id VARCHAR, ko_raw VARCHAR",
        "SELECT DISTINCT i.organism_code, i.source_path, trim(r.gene_id) AS gene_id, trim(r.ko_raw) AS ko_raw",
    ),
    "gene_pathway": (
        "_kegg_gene_pathway",
        "gene_id VARCHAR, pathway_raw VARCHAR",
        "SELECT DISTINCT i.organism_code, i.source_path, trim(r.gene_id) AS gene_id, trim(r.pathway_raw) AS pathway_raw",
    ),
}


def write_native_mapping_publication(
    snapshot: MappingSnapshot,
    path: Path,
    *,
    if_exists: Literal["fail", "replace"],
) -> DuckDBWriteResult:
    """Build one mapping publication with native DuckDB scans and SQL joins."""
    destination = preflight_publication_destination(path, if_exists=if_exists)
    package_version = require_package_version()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = _create_stage_path(destination, suffix=".duckdb")
    stage.unlink()
    connection: duckdb.DuckDBPyConnection | None = None
    spill_directory: Any = None
    row_counts = dict.fromkeys(TABLE_SCHEMAS, 0)
    sources: list[SourceFileRecord] = []
    try:
        spill_directory = tempfile.TemporaryDirectory(
            prefix="bioextract-kegg-duckdb-",
            dir=str(destination.parent),
        )
        memory_limit, max_temp_directory_size = _resource_limits(
            Path(spill_directory.name)
        )
        connection = _connect_publication(
            stage,
            threads=min(pl.thread_pool_size(), _THREAD_CAP),
            temp_directory=Path(spill_directory.name),
            memory_limit=memory_limit,
            max_temp_directory_size=max_temp_directory_size,
        )
        _create_metadata_schema(connection)
        _create_biological_tables(connection)

        work = resolve_organism_work(
            snapshot,
            validate_role_files=snapshot.mode == "directory",
        )
        inventory, role_paths = _build_inventory(snapshot, work)
        _create_input_file(connection, inventory)
        sources.extend(
            SourceFileRecord(
                logical_name=logical_name,
                path=source_path,
                media_type=MEDIA_TYPE_TSV,
                bytes=None,
                sha256=None,
            )
            for logical_name, _, _, source_path in inventory
        )

        for role in ("organism_list", "ko_pathway"):
            paths = role_paths.get(role, ())
            if not paths:
                continue
            _scan_role(connection, role, paths)
            _validate_role_content(connection, {role: paths})
            if role == "organism_list":
                _create_organism_meta(connection)
            else:
                _create_ko_pathway(connection)
            connection.execute(f'DROP TABLE IF EXISTS "_kegg_raw_{role}"')
        capabilities = source_capabilities(snapshot)
        _build_organism_table(connection, capabilities)
        _build_gene_table(connection, capabilities)
        _build_ko_table(connection, capabilities)

        row_counts = {
            table_name: _scalar_int(connection, f'SELECT count(*) FROM "{table_name}"')
            for table_name in TABLE_SCHEMAS
        }
        metadata = _publication_metadata(
            package_version=package_version,
            resource_name="kegg",
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile=SOURCE_SCHEMA_PROFILE,
            source_schema_version=None,
            scope="mapping",
            release_version=snapshot.release_version,
            release_version_source=(
                "caller" if snapshot.release_version is not None else None
            ),
            validation_issue_count=0,
        )
        metadata.update(
            {
                f"bioextract.capability.{name}": str(value).lower()
                for name, value in capabilities.items()
            }
        )
        metadata["bioextract.organism_scope_mode"] = (
            "selected"
            if snapshot.organism_scope is not None or snapshot.mode == "files"
            else "all_available_at_build"
        )
        relations = _relation_specs()
        _write_duckdb_metadata(
            connection,
            metadata=metadata,
            sources=sources,
            relations=relations,
            row_counts=row_counts,
            column_mappings=(),
            validation_issues=(),
        )
        _drop_private_tables(connection)
        connection.execute("CHECKPOINT")
        connection.close()
        connection = None
        _validate_duckdb_publication(stage, relations=relations, row_counts=row_counts)
        # Mapping-specific validation is imported lazily to avoid a module cycle.
        from .publication import validate_mapping_publication

        validate_mapping_publication(stage)
        stage.replace(destination)
    except BaseException:
        if connection is not None:
            connection.close()
        raise
    finally:
        if spill_directory is not None:
            spill_directory.cleanup()
        stage.unlink(missing_ok=True)
    return DuckDBWriteResult(
        path=destination,
        resource_name="kegg",
        resource_schema_version=SCHEMA_VERSION,
        tables=tuple(TABLE_SCHEMAS),
        row_counts=row_counts,
        validation_issue_count=0,
    )


def _build_inventory(
    snapshot: MappingSnapshot,
    work: Sequence[tuple[str, Mapping[str, Path]]],
) -> tuple[
    list[tuple[str, str | None, str, Path]],
    dict[str, list[Path]],
]:
    inventory: list[tuple[str, str | None, str, Path]] = []
    role_paths: dict[str, list[Path]] = {}
    seen_paths: set[Path] = set()
    for organism_code, roles in work:
        for role, path in roles.items():
            resolved = path.resolve()
            if resolved in seen_paths:
                continue
            seen_paths.add(resolved)
            inventory.append(
                (f"organism/{organism_code}/{role}", organism_code, role, resolved)
            )
            role_paths.setdefault(role, []).append(resolved)
    globals_by_role = {
        "organism_list": snapshot.organism_list,
        "ko_pathway": snapshot.ko_pathway,
    }
    for role, value in globals_by_role.items():
        if value is None:
            continue
        resolved = value.resolve()
        if resolved in seen_paths:
            continue
        seen_paths.add(resolved)
        inventory.append((f"global/{role}", None, role, resolved))
        role_paths.setdefault(role, []).append(resolved)
    return inventory, role_paths


def _create_input_file(
    connection: duckdb.DuckDBPyConnection,
    inventory: Sequence[tuple[str, str | None, str, Path]],
) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_input_file (
            logical_name VARCHAR,
            organism_code VARCHAR,
            role VARCHAR,
            source_path VARCHAR
        )
        """
    )
    connection.executemany(
        "INSERT INTO _kegg_input_file VALUES (?, ?, ?, ?)",
        [
            (logical_name, organism_code, role, str(path))
            for logical_name, organism_code, role, path in inventory
        ],
    )


def _scan_role(
    connection: duckdb.DuckDBPyConnection,
    role: str,
    paths: Sequence[Path],
) -> None:
    columns = _ROLE_COLUMNS[role]
    column_map = "{" + ", ".join(f"'{name}': 'VARCHAR'" for name in columns) + "}"
    raw_table = f"_kegg_raw_{role}"
    rejects_table = f"_kegg_reject_{role}_errors"
    rejects_scan = f"_kegg_reject_{role}_scan"
    path_values = [str(path) for path in paths]
    query = f'''
        CREATE TEMP TABLE "{raw_table}" AS
        SELECT * FROM read_csv(
            ?,
            auto_detect=false,
            header=false,
            delim='\\t',
            quote='',
            escape='',
            columns={column_map},
            filename=true,
            union_by_name=false,
            null_padding=false,
            strict_mode=true,
                store_rejects=true,
                rejects_table='{rejects_table}',
                rejects_scan='{rejects_scan}',
                rejects_limit={_REJECTS_LIMIT}
        )
    '''
    try:
        connection.execute(query, [path_values])
    except duckdb.Error as error:
        raise IntegrityError(
            f"KEGG mapping parse error: role={role!r}, "
            f"path={paths[0]}, cause=engine_parse"
        ) from error
    try:
        rejected = connection.execute(
            f'SELECT count(*) FROM "{rejects_table}"'
        ).fetchone()
    except duckdb.Error:
        rejected = None
    if rejected is not None and int(rejected[0]) > 0:
        raise IntegrityError(
            f"KEGG mapping parse error: role={role!r}, "
            f"path={paths[0]}, cause=csv_rejects:{int(rejected[0])}"
        )


def _validate_role_content(
    connection: duckdb.DuckDBPyConnection,
    role_paths: Mapping[str, Sequence[Path]],
) -> None:
    for role, paths in role_paths.items():
        if not paths:
            continue
        if role == "gene_list":
            _assert_valid(
                connection,
                "_kegg_raw_gene_list",
                "gene_list",
                "gene_id = ''",
            )
        elif role in {"uniprot_conversion", "ncbi_gene_conversion"}:
            _assert_valid(
                connection,
                f"_kegg_raw_{role}",
                role,
                "xref = '' OR gene_id = ''",
            )
        elif role == "gene_ko":
            _assert_valid(
                connection,
                "_kegg_raw_gene_ko",
                role,
                "gene_id = '' OR ko_raw = ''",
            )
        elif role == "gene_pathway":
            _assert_valid(
                connection,
                "_kegg_raw_gene_pathway",
                role,
                "gene_id = '' OR pathway_raw = ''",
            )
    if role_paths.get("organism_list"):
        _assert_valid(
            connection,
            "_kegg_raw_organism_list",
            "global/organism_list",
            "genome_id = '' OR organism_code = ''",
        )
    if role_paths.get("ko_pathway"):
        _assert_valid(
            connection,
            "_kegg_raw_ko_pathway",
            "global/ko_pathway",
            "ko_raw = '' OR pathway_raw = ''",
        )


def _assert_valid(
    connection: duckdb.DuckDBPyConnection,
    table: str,
    role: str,
    predicate: str,
    *,
    path_column: str = "filename",
) -> None:
    row = connection.execute(
        f'SELECT "{path_column}" FROM "{table}" WHERE {predicate} LIMIT 1'
    ).fetchone()
    if row is not None:
        raise IntegrityError(
            f"KEGG mapping parse error: role={role!r}, "
            f"path={row[0]}, cause=invalid_required_field"
        )


def _build_organism_table(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
) -> None:
    if capabilities.get("organism_list"):
        metadata_join = "LEFT JOIN _kegg_organism_meta AS m USING (organism_code)"
        genome_id = "m.genome_id"
        organism_name = "m.organism_name"
        taxonomy = "coalesce(m.taxonomy_lineage, []::VARCHAR[])"
    else:
        metadata_join = ""
        genome_id = "NULL::VARCHAR"
        organism_name = "NULL::VARCHAR"
        taxonomy = "NULL::VARCHAR[]"
    connection.execute(
        f"""
        INSERT INTO organism
        SELECT DISTINCT
            i.organism_code,
            {genome_id},
            {organism_name},
            {taxonomy}
        FROM _kegg_input_file AS i
        {metadata_join}
        WHERE i.organism_code IS NOT NULL
        ORDER BY i.organism_code
        """
    )


def _create_organism_meta(connection: duckdb.DuckDBPyConnection) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_organism_meta AS
        SELECT DISTINCT
            trim(organism_code) AS organism_code,
            trim(genome_id) AS genome_id,
            nullif(trim(organism_name), '') AS organism_name,
            list_transform(
                list_filter(string_split(taxonomy, ';'), x -> trim(x) <> ''),
                x -> trim(x)
            ) AS taxonomy_lineage
        FROM _kegg_raw_organism_list
        """
    )
    _assert_valid(
        connection,
        "_kegg_raw_organism_list",
        "global/organism_list",
        "NOT regexp_matches(genome_id, '^T[0-9]+$')",
    )
    conflict = connection.execute(
        """
        SELECT organism_code
        FROM _kegg_organism_meta
        GROUP BY organism_code
        HAVING count(*) > 1
        LIMIT 1
        """
    ).fetchone()
    if conflict is not None:
        raise IntegrityError(
            "KEGG mapping parse error: role='global/organism_list', "
            f"cause=conflicting_metadata:{conflict[0]}"
        )


def _build_gene_table(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_ko_ids (ko_id VARCHAR)
        """
    )
    for partition in range(_AGGREGATION_PARTITIONS):
        role_paths = _partition_role_paths(connection, partition)
        for role in _ORGANISM_ROLES:
            paths = role_paths.get(role, ())
            if not paths:
                continue
            _scan_role(connection, role, paths)
            _validate_role_content(connection, {role: paths})
            _create_normalized_gene_role(connection, role)
            connection.execute(f'DROP TABLE IF EXISTS "_kegg_raw_{role}"')
        _create_normalized_gene_tables(connection)
        _validate_normalized_gene_tables(connection)
        _create_gene_attributes(connection, capabilities)
        connection.execute(
            """
            CREATE TEMP TABLE _kegg_gene_universe AS
            SELECT DISTINCT organism_code, gene_id AS kegg_gene_id
            FROM (
                SELECT organism_code, gene_id FROM _kegg_gene_list
                UNION ALL SELECT organism_code, gene_id FROM _kegg_uniprot
                UNION ALL SELECT organism_code, gene_id FROM _kegg_ncbi
                UNION ALL SELECT organism_code, gene_id FROM _kegg_gene_ko
                UNION ALL SELECT organism_code, gene_id FROM _kegg_gene_pathway
            )
            """
        )
        if capabilities.get("gene_ko"):
            connection.execute(
                """
                INSERT INTO _kegg_ko_ids
                SELECT DISTINCT substr(ko_raw, 4) FROM _kegg_gene_ko
                """
            )
        _create_gene_aggregates(connection, capabilities, partition=partition)
        _insert_gene_partition(connection, capabilities, partition=partition)
        for table in (
            "_kegg_gene_list",
            "_kegg_uniprot",
            "_kegg_ncbi",
            "_kegg_gene_ko",
            "_kegg_gene_pathway",
            "_kegg_gene_attributes",
            "_kegg_gene_universe",
            "_kegg_uniprot_agg",
            "_kegg_ncbi_agg",
            "_kegg_ko_agg",
            "_kegg_pathway_agg",
        ):
            connection.execute(f'DROP TABLE IF EXISTS "{table}"')


def _insert_gene_partition(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
    *,
    partition: int,
) -> None:
    empty_aliases = "[]::VARCHAR[]"
    empty_uniprot = "[]::STRUCT(uniprot_id VARCHAR)[]"
    empty_ncbi = "[]::STRUCT(ncbi_gene_id VARCHAR)[]"
    empty_ko = "[]::STRUCT(ko_id VARCHAR)[]"
    empty_pathway = "[]::STRUCT(kegg_pathway_id VARCHAR, pathway_map_id VARCHAR)[]"
    predicate = f"hash(u.organism_code) % {_AGGREGATION_PARTITIONS} = {partition}"
    connection.execute(
        f"""
        INSERT INTO gene_annotation
        SELECT
            u.organism_code,
            u.kegg_gene_id,
            a.gene_type,
            a.genomic_position,
            a.gene_symbol,
            {f"coalesce(a.gene_aliases, {empty_aliases})" if capabilities.get("gene_list") else "NULL::VARCHAR[]"},
            a.gene_description,
            {f"coalesce(up.mappings, {empty_uniprot})" if capabilities.get("uniprot_conversion") else "NULL::STRUCT(uniprot_id VARCHAR)[]"},
            {f"coalesce(nc.mappings, {empty_ncbi})" if capabilities.get("ncbi_gene_conversion") else "NULL::STRUCT(ncbi_gene_id VARCHAR)[]"},
            {f"coalesce(ko.mappings, {empty_ko})" if capabilities.get("gene_ko") else "NULL::STRUCT(ko_id VARCHAR)[]"},
            {f"coalesce(pw.mappings, {empty_pathway})" if capabilities.get("gene_pathway") else "NULL::STRUCT(kegg_pathway_id VARCHAR, pathway_map_id VARCHAR)[]"}
        FROM _kegg_gene_universe AS u
        LEFT JOIN _kegg_gene_attributes AS a
          ON a.organism_code = u.organism_code AND a.gene_id = u.kegg_gene_id
        LEFT JOIN _kegg_uniprot_agg AS up
          ON up.organism_code = u.organism_code AND up.gene_id = u.kegg_gene_id
        LEFT JOIN _kegg_ncbi_agg AS nc
          ON nc.organism_code = u.organism_code AND nc.gene_id = u.kegg_gene_id
        LEFT JOIN _kegg_ko_agg AS ko
          ON ko.organism_code = u.organism_code AND ko.gene_id = u.kegg_gene_id
        LEFT JOIN _kegg_pathway_agg AS pw
          ON pw.organism_code = u.organism_code AND pw.gene_id = u.kegg_gene_id
        WHERE {predicate}
        ORDER BY u.organism_code, u.kegg_gene_id
        """
    )


def _partition_role_paths(
    connection: duckdb.DuckDBPyConnection,
    partition: int,
) -> dict[str, tuple[Path, ...]]:
    rows = connection.execute(
        f"""
        SELECT role, source_path
        FROM _kegg_input_file
        WHERE organism_code IS NOT NULL
          AND hash(organism_code) % {_AGGREGATION_PARTITIONS} = {partition}
        ORDER BY role, source_path
        """
    ).fetchall()
    paths: dict[str, list[Path]] = {}
    for role, source_path in rows:
        paths.setdefault(str(role), []).append(Path(str(source_path)))
    return {role: tuple(values) for role, values in paths.items()}


def _create_normalized_gene_tables(connection: duckdb.DuckDBPyConnection) -> None:
    for role, (table, fields, _select) in _NORMALIZED_GENE_DEFINITIONS.items():
        if _temp_table_exists(connection, table):
            continue
        raw_table = f"_kegg_raw_{role}"
        if _temp_table_exists(connection, raw_table):
            _create_normalized_gene_role(connection, role)
        else:
            connection.execute(
                f"""CREATE TEMP TABLE \"{table}\" (
                    organism_code VARCHAR, source_path VARCHAR, {fields}
                )"""
            )


def _create_normalized_gene_role(
    connection: duckdb.DuckDBPyConnection,
    role: str,
) -> None:
    table, _fields, select = _NORMALIZED_GENE_DEFINITIONS[role]
    raw_table = f"_kegg_raw_{role}"
    connection.execute(
        f"""
        CREATE TEMP TABLE \"{table}\" AS
        {select}
        FROM \"{raw_table}\" AS r
        JOIN _kegg_input_file AS i
          ON i.source_path = r.filename AND i.role = '{role}'
        """
    )


def _validate_normalized_gene_tables(connection: duckdb.DuckDBPyConnection) -> None:
    checks = (
        (
            "_kegg_gene_list",
            "gene_list",
            "gene_id = '' OR NOT starts_with(gene_id, organism_code || ':')",
        ),
        (
            "_kegg_uniprot",
            "uniprot_conversion",
            "xref = '' OR strpos(xref, ':') <= 1 OR gene_id = '' OR NOT starts_with(gene_id, organism_code || ':')",
        ),
        (
            "_kegg_ncbi",
            "ncbi_gene_conversion",
            "NOT regexp_matches(xref, '^ncbi-geneid:[0-9]+$') OR gene_id = '' OR NOT starts_with(gene_id, organism_code || ':')",
        ),
        (
            "_kegg_gene_ko",
            "gene_ko",
            "gene_id = '' OR NOT starts_with(gene_id, organism_code || ':') OR NOT regexp_matches(ko_raw, '^ko:K[0-9]{5}$')",
        ),
        (
            "_kegg_gene_pathway",
            "gene_pathway",
            "gene_id = '' OR NOT starts_with(gene_id, organism_code || ':') OR NOT regexp_matches(pathway_raw, '^path:' || organism_code || '[0-9]{5}$')",
        ),
    )
    for table, role, predicate in checks:
        _assert_valid(connection, table, role, predicate, path_column="source_path")
    if _scalar_int(
        connection,
        """
        SELECT count(*) FROM (
            SELECT organism_code, gene_id
            FROM _kegg_gene_list
            GROUP BY organism_code, gene_id
            HAVING count(DISTINCT (gene_type, genomic_position, display)) > 1
        )
        """,
    ):
        raise IntegrityError(
            "KEGG mapping parse error: role='gene_list', cause=conflicting_metadata"
        )


def _create_gene_attributes(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
) -> None:
    if not capabilities.get("gene_list"):
        connection.execute(
            """
            CREATE TEMP TABLE _kegg_gene_attributes AS
            SELECT NULL::VARCHAR AS organism_code, NULL::VARCHAR AS gene_id,
                NULL::VARCHAR AS gene_type, NULL::VARCHAR AS genomic_position,
                NULL::VARCHAR AS gene_symbol, NULL::VARCHAR[] AS gene_aliases,
                NULL::VARCHAR AS gene_description
            WHERE false
            """
        )
        return
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_gene_attributes AS
        WITH named AS (
            SELECT DISTINCT organism_code, gene_id, gene_type, genomic_position,
                display,
                list_transform(
                    list_filter(
                        string_split(
                            CASE WHEN strpos(display, ';') = 0
                                 THEN display
                                 ELSE substr(display, 1, strpos(display, ';') - 1)
                            END,
                            ','
                        ),
                        x -> trim(x) <> ''
                    ),
                    x -> trim(x)
                ) AS name_parts
            FROM _kegg_gene_list
        )
        SELECT organism_code, gene_id, gene_type, genomic_position,
            CASE WHEN strpos(display, ';') = 0 THEN NULL
                 ELSE nullif(list_extract(name_parts, 1), '') END AS gene_symbol,
            CASE WHEN strpos(display, ';') = 0 THEN []::VARCHAR[]
                 ELSE list_sort(list_distinct(list_slice(
                     name_parts, 2, length(name_parts)
                 ))) END AS gene_aliases,
            CASE WHEN strpos(display, ';') = 0
                 THEN nullif(trim(display), '')
                 ELSE nullif(trim(substr(display, strpos(display, ';') + 1)), '')
            END AS gene_description
        FROM named
        """
    )


def _create_gene_aggregates(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
    *,
    partition: int,
) -> None:
    partition_predicate = (
        f"hash(organism_code) % {_AGGREGATION_PARTITIONS} = {partition}"
    )
    if capabilities.get("uniprot_conversion"):
        connection.execute(
            f"""
            CREATE TEMP TABLE _kegg_uniprot_agg AS
            SELECT organism_code, gene_id,
                list(struct_pack(uniprot_id := accession) ORDER BY accession) AS mappings
            FROM (
                SELECT DISTINCT organism_code, gene_id,
                    substr(xref, strpos(xref, ':') + 1) AS accession
                FROM _kegg_uniprot
                WHERE {partition_predicate}
            )
            GROUP BY organism_code, gene_id
            """
        )
    else:
        _empty_aggregate(connection, "_kegg_uniprot_agg", "uniprot_id VARCHAR")
    if capabilities.get("ncbi_gene_conversion"):
        connection.execute(
            f"""
            CREATE TEMP TABLE _kegg_ncbi_agg AS
            SELECT organism_code, gene_id,
                list(struct_pack(ncbi_gene_id := accession) ORDER BY accession) AS mappings
            FROM (
                SELECT DISTINCT organism_code, gene_id,
                    substr(xref, strpos(xref, ':') + 1) AS accession
                FROM _kegg_ncbi
                WHERE {partition_predicate}
            )
            GROUP BY organism_code, gene_id
            """
        )
    else:
        _empty_aggregate(connection, "_kegg_ncbi_agg", "ncbi_gene_id VARCHAR")
    if capabilities.get("gene_ko"):
        connection.execute(
            f"""
            CREATE TEMP TABLE _kegg_ko_agg AS
            SELECT organism_code, gene_id,
                list(struct_pack(ko_id := ko_id) ORDER BY ko_id) AS mappings
            FROM (
                SELECT DISTINCT organism_code, gene_id, substr(ko_raw, 4) AS ko_id
                FROM _kegg_gene_ko
                WHERE {partition_predicate}
            )
            GROUP BY organism_code, gene_id
            """
        )
    else:
        _empty_aggregate(connection, "_kegg_ko_agg", "ko_id VARCHAR")
    if capabilities.get("gene_pathway"):
        connection.execute(
            f"""
            CREATE TEMP TABLE _kegg_pathway_agg AS
            SELECT organism_code, gene_id,
                list(struct_pack(
                    kegg_pathway_id := pathway_id,
                    pathway_map_id := concat('map', right(pathway_id, 5))
                ) ORDER BY pathway_id) AS mappings
            FROM (
                SELECT DISTINCT organism_code, gene_id,
                    substr(pathway_raw, 6) AS pathway_id
                FROM _kegg_gene_pathway
                WHERE {partition_predicate}
            )
            GROUP BY organism_code, gene_id
            """
        )
    else:
        _empty_aggregate(
            connection,
            "_kegg_pathway_agg",
            "kegg_pathway_id VARCHAR, pathway_map_id VARCHAR",
        )


def _empty_aggregate(
    connection: duckdb.DuckDBPyConnection,
    table: str,
    struct_fields: str,
    *,
    key_fields: str = "organism_code VARCHAR, gene_id VARCHAR",
) -> None:
    connection.execute(
        f"""
        CREATE TEMP TABLE \"{table}\" (
            {key_fields},
            mappings STRUCT({struct_fields})[]
        )
        """
    )


def _build_ko_table(
    connection: duckdb.DuckDBPyConnection,
    capabilities: Mapping[str, bool],
) -> None:
    if capabilities.get("ko_pathway"):
        if not _temp_table_exists(connection, "_kegg_ko_pathway"):
            _create_empty_ko_pathway(connection)
        _assert_valid(
            connection,
            "_kegg_ko_pathway",
            "global/ko_pathway",
            "NOT regexp_matches(ko_raw, '^ko:K[0-9]{5}$') OR NOT regexp_matches(pathway_raw, '^path:(map|ko)[0-9]{5}$')",
            path_column="source_path",
        )
    else:
        connection.execute(
            """
            CREATE TEMP TABLE _kegg_ko_pathway (
                source_path VARCHAR,
                ko_raw VARCHAR,
                pathway_raw VARCHAR
            )
            """
        )
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_ko_universe AS
        SELECT DISTINCT ko_id
        FROM (
            SELECT ko_id FROM _kegg_ko_ids
            UNION ALL
            SELECT substr(ko_raw, 4) AS ko_id FROM _kegg_ko_pathway
        )
        WHERE ko_id IS NOT NULL AND ko_id <> ''
        """
    )
    for partition in range(_AGGREGATION_PARTITIONS):
        if capabilities.get("ko_pathway"):
            connection.execute(
                f"""
                CREATE TEMP TABLE _kegg_ko_pathway_agg AS
                SELECT ko_id,
                    list(struct_pack(
                        kegg_pathway_id := pathway_id,
                        pathway_namespace := namespace,
                        pathway_map_id := concat('map', right(pathway_id, 5))
                    ) ORDER BY pathway_id) AS mappings
                FROM (
                    SELECT DISTINCT substr(ko_raw, 4) AS ko_id,
                        substr(pathway_raw, 6) AS pathway_id,
                        CASE
                            WHEN starts_with(pathway_raw, 'path:ko') THEN 'ko'
                            ELSE 'map'
                        END AS namespace
                    FROM _kegg_ko_pathway
                    WHERE hash(substr(ko_raw, 4)) % {_AGGREGATION_PARTITIONS} = {partition}
                )
                GROUP BY ko_id
                """
            )
        else:
            _empty_aggregate(
                connection,
                "_kegg_ko_pathway_agg",
                "kegg_pathway_id VARCHAR, pathway_namespace VARCHAR, pathway_map_id VARCHAR",
                key_fields="ko_id VARCHAR",
            )
        empty_pathways = "[]::STRUCT(kegg_pathway_id VARCHAR, pathway_namespace VARCHAR, pathway_map_id VARCHAR)[]"
        pathway_expr = (
            f"coalesce(p.mappings, {empty_pathways})"
            if capabilities.get("ko_pathway")
            else "NULL::STRUCT(kegg_pathway_id VARCHAR, pathway_namespace VARCHAR, pathway_map_id VARCHAR)[]"
        )
        connection.execute(
            f"""
            INSERT INTO ko_annotation
            SELECT u.ko_id, {pathway_expr}
            FROM _kegg_ko_universe AS u
            LEFT JOIN _kegg_ko_pathway_agg AS p USING (ko_id)
            WHERE hash(u.ko_id) % {_AGGREGATION_PARTITIONS} = {partition}
            ORDER BY u.ko_id
            """
        )
        connection.execute('DROP TABLE IF EXISTS "_kegg_ko_pathway_agg"')


def _create_ko_pathway(connection: duckdb.DuckDBPyConnection) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_ko_pathway AS
        SELECT DISTINCT i.source_path, trim(r.ko_raw) AS ko_raw,
            trim(r.pathway_raw) AS pathway_raw
        FROM _kegg_raw_ko_pathway AS r
        JOIN _kegg_input_file AS i
          ON i.source_path = r.filename AND i.role = 'ko_pathway'
        """
    )


def _create_empty_ko_pathway(connection: duckdb.DuckDBPyConnection) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _kegg_ko_pathway (
            source_path VARCHAR,
            ko_raw VARCHAR,
            pathway_raw VARCHAR
        )
        """
    )


def _drop_private_tables(connection: duckdb.DuckDBPyConnection) -> None:
    names = connection.execute(
        """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = 'main'
          AND table_name LIKE '_kegg_%'
        """
    ).fetchall()
    for (name,) in names:
        connection.execute(f'DROP TABLE IF EXISTS "{name}"')
    names = connection.execute(
        """
        SELECT table_schema, table_name
        FROM information_schema.tables
        WHERE table_name LIKE '_kegg_reject_%'
        """
    ).fetchall()
    for schema, name in names:
        connection.execute(f'DROP TABLE IF EXISTS "{schema}"."{name}"')


def _temp_table_exists(connection: duckdb.DuckDBPyConnection, name: str) -> bool:
    row = connection.execute(
        "SELECT 1 FROM information_schema.tables "
        "WHERE table_schema = 'main' AND table_name = ?",
        [name],
    ).fetchone()
    return row is not None


def _resource_limits(temp_directory: Path) -> tuple[str, str]:
    """Return private memory and spill bounds for one staged build."""
    memory_bytes = max(1 << 30, _available_memory_bytes() // 2)
    free_bytes = shutil.disk_usage(temp_directory).free
    temp_bytes = max(64 << 20, free_bytes // 2)
    return f"{memory_bytes}B", f"{temp_bytes}B"


def _available_memory_bytes() -> int:
    """Prefer the cgroup limit, falling back to host physical memory."""
    try:
        host_bytes = int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (AttributeError, OSError, ValueError):
        host_bytes = 0
    for candidate in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            value = candidate.read_text(encoding="ascii").strip()
        except OSError:
            continue
        if value.isdigit() and int(value) > 0:
            limit = int(value)
            if host_bytes <= 0 or limit <= host_bytes:
                return limit
    return host_bytes or (8 << 30)


def _create_biological_tables(connection: duckdb.DuckDBPyConnection) -> None:
    for table_name, schema in TABLE_SCHEMAS.items():
        columns = ", ".join(
            f'"{name}" {_duckdb_type(dtype)}' for name, dtype in schema.items()
        )
        connection.execute(f'CREATE TABLE "{table_name}" ({columns})')


def _relation_specs() -> tuple[RelationSpec, ...]:
    return tuple(
        RelationSpec(
            table_name=table_name,
            frame=_empty_lazy(schema),
            role=_TABLE_ROLES[table_name],
        )
        for table_name, schema in TABLE_SCHEMAS.items()
    )


def _empty_lazy(schema: Mapping[str, Any]) -> Any:
    import polars as pl

    return pl.DataFrame(schema=schema).lazy()


def _duckdb_type(dtype: Any) -> str:
    import polars as pl

    if dtype == pl.String:
        return "VARCHAR"
    if isinstance(dtype, pl.List):
        return f"{_duckdb_type(dtype.inner)}[]"
    if isinstance(dtype, pl.Struct):
        fields = ", ".join(
            f"{field.name} {_duckdb_type(field.dtype)}" for field in dtype.fields
        )
        return f"STRUCT({fields})"
    raise TypeError(f"Unsupported KEGG mapping type: {dtype!r}")


def _scalar_int(connection: duckdb.DuckDBPyConnection, query: str) -> int:
    row = connection.execute(query).fetchone()
    if row is None:
        raise IntegrityError(f"KEGG mapping validation query returned no row: {query}")
    return int(row[0])
