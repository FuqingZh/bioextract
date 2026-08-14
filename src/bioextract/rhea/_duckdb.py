from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable, Mapping
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Literal

import duckdb
import polars as pl

from bioextract._publication import (
    METADATA_SCHEMA_VERSION,
    validate_duckdb_metadata_v2,
)

from .constant import SCHEMA_VERSION, SOURCE_SCHEMA_PROFILE
from .rhea import RheaWriteResult
from .util import (
    coefficient_numeric,
    compound_type,
    infer_media_type,
    is_gzip_file,
    iter_sdf_structures,
    parse_rdf_snapshot,
    public_compound_accession,
    reaction_direction,
    read_release_properties,
)


def write_rhea_duckdb(
    *,
    sources: Mapping[str, Path],
    display_paths: Mapping[str, str],
    scope: str,
    path: Path,
    if_exists: Literal["fail", "replace"],
) -> RheaWriteResult:
    """Build a Rhea database in staging and atomically commit it."""

    if path.exists() and if_exists == "fail":
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    descriptor, staging_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(descriptor)
    staging = Path(staging_name)
    staging.unlink()

    release_number, release_date = _release_metadata(sources)
    connection: duckdb.DuckDBPyConnection | None = None
    try:
        connection = _connect_publication(staging)
        _create_metadata_tables(connection)

        if "rdf" in sources:
            _load_rdf(connection, sources["rdf"])
        if "directions" in sources:
            _load_directions(connection, sources["directions"])
            if "rdf" in sources:
                _validate_rdf_directions(connection)
        if "relationships" in sources:
            _load_relationships(connection, sources["relationships"])
        if "obsoletes" in sources:
            _load_obsoletes(connection, sources["obsoletes"])
        if "reaction_smiles" in sources:
            _load_reaction_smiles(connection, sources["reaction_smiles"])

        if "sdf" in sources:
            _load_sdf(connection, sources["sdf"])
        if "chebi_names" in sources:
            _load_chebi_names(connection, sources["chebi_names"])
        if "chebi_ph7_3_mapping" in sources:
            _load_chebi_mapping(connection, sources["chebi_ph7_3_mapping"])

        if "xrefs" in sources:
            _load_xrefs(connection, sources["xrefs"])
        if "uniprot_sprot" in sources or "uniprot_trembl" in sources:
            _load_uniprot(
                connection,
                file_sprot=sources.get("uniprot_sprot"),
                file_trembl=sources.get("uniprot_trembl"),
            )

        if "xrefs" in sources and "ec" in sources:
            _validate_specialized_xref(connection, sources["ec"], database_name="EC")
        if "xrefs" in sources and "go" in sources:
            _validate_specialized_xref(connection, sources["go"], database_name="GO")

        _create_views(connection)
        _record_metadata(
            connection,
            sources=sources,
            display_paths=display_paths,
            scope=scope,
            release_number=release_number,
            release_date=release_date,
        )
        row_counts = _record_table_inventory(connection)
        connection.close()
        connection = None

        _validate_staged_database(staging, row_counts=row_counts)
        os.replace(staging, path)
        return RheaWriteResult(
            path=path,
            scope=scope,
            tables=tuple(sorted(row_counts)),
            row_counts=row_counts,
            source_files=dict(display_paths),
            release_number=release_number,
            release_date=release_date,
        )
    except BaseException:
        if connection is not None:
            connection.close()
        _remove_staging_files(staging)
        raise


def _create_metadata_tables(connection: duckdb.DuckDBPyConnection) -> None:
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
            bytes UBIGINT,
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


def _load_rdf(connection: duckdb.DuckDBPyConnection, file_rdf: Path) -> None:
    snapshot = parse_rdf_snapshot(file_rdf)
    connection.execute(
        """
        CREATE TABLE reaction (
            rhea_id BIGINT PRIMARY KEY,
            master_id BIGINT,
            direction VARCHAR,
            accession VARCHAR,
            equation VARCHAR,
            equation_html VARCHAR,
            status VARCHAR,
            is_balanced BOOLEAN,
            is_transport BOOLEAN,
            comment VARCHAR,
            is_obsolete BOOLEAN NOT NULL
        )
        """
    )
    _insert_rows_batched(
        connection,
        "reaction",
        [
            (
                reaction.rhea_id,
                *reaction_direction(reaction),
                reaction.accession,
                reaction.equation,
                reaction.equation_html,
                reaction.status,
                reaction.is_balanced,
                reaction.is_transport,
                reaction.comment,
                (reaction.status or "").lower() == "obsolete",
            )
            for reaction in snapshot.reactions.values()
        ],
    )
    connection.execute(
        """
        CREATE TABLE reaction_side (
            side_id VARCHAR PRIMARY KEY,
            master_id BIGINT NOT NULL,
            side VARCHAR NOT NULL,
            curated_order INTEGER
        )
        """
    )
    _insert_rows_batched(
        connection,
        "reaction_side",
        [
            (side.side_id, side.master_id, side.side, side.curated_order)
            for side in snapshot.sides.values()
        ],
    )
    connection.execute(
        """
        CREATE TABLE compound (
            compound_id VARCHAR PRIMARY KEY,
            rhea_compound_id BIGINT,
            source_accession VARCHAR,
            public_accession VARCHAR,
            compound_type VARCHAR,
            name VARCHAR,
            name_html VARCHAR,
            formula VARCHAR,
            charge_text VARCHAR,
            charge_numeric INTEGER,
            chebi_id VARCHAR,
            underlying_chebi_id VARCHAR,
            polymerization_index VARCHAR
        )
        """
    )
    _insert_rows_batched(
        connection,
        "compound",
        [
            (
                compound.compound_id,
                compound.rhea_compound_id,
                compound.source_accession,
                public_compound_accession(compound.source_accession),
                compound_type(compound),
                compound.name,
                compound.name_html,
                compound.formula,
                compound.charge,
                _optional_int(compound.charge),
                _chebi_curie(compound.chebi_id),
                _chebi_curie(compound.underlying_chebi_id),
                compound.polymerization_index,
            )
            for compound in snapshot.compounds.values()
        ],
    )
    connection.execute(
        """
        CREATE TABLE compound_reactive_part (
            compound_id VARCHAR NOT NULL,
            reactive_part_id VARCHAR NOT NULL,
            position VARCHAR,
            PRIMARY KEY (compound_id, reactive_part_id)
        )
        """
    )
    _insert_rows_batched(
        connection,
        "compound_reactive_part",
        [
            (
                compound.compound_id,
                reactive_part_id,
                snapshot.compounds[reactive_part_id].position
                if reactive_part_id in snapshot.compounds
                else None,
            )
            for compound in snapshot.compounds.values()
            for reactive_part_id in sorted(compound.reactive_part_ids)
        ],
    )
    connection.execute(
        """
        CREATE TABLE reaction_participant (
            side_id VARCHAR NOT NULL,
            participant_id VARCHAR NOT NULL,
            master_id BIGINT NOT NULL,
            side VARCHAR NOT NULL,
            compound_id VARCHAR,
            coefficient_text VARCHAR NOT NULL,
            coefficient_numeric DOUBLE,
            location VARCHAR,
            PRIMARY KEY (side_id, participant_id)
        )
        """
    )
    participant_rows: list[
        tuple[str, str, int, str, str | None, str, float | None, str | None]
    ] = []
    for membership in snapshot.memberships:
        side = snapshot.sides[membership.side_id]
        participant = snapshot.participants.get(membership.participant_id)
        coefficient = snapshot.coefficient_by_predicate.get(
            membership.coefficient_predicate, "1"
        )
        participant_rows.append(
            (
                membership.side_id,
                membership.participant_id,
                side.master_id,
                side.side,
                None if participant is None else participant.compound_id,
                coefficient,
                coefficient_numeric(coefficient),
                None if participant is None else participant.location,
            )
        )
    _insert_rows_batched(connection, "reaction_participant", participant_rows)
    connection.execute(
        """
        CREATE TABLE reaction_publication (
            rhea_id BIGINT NOT NULL,
            pubmed_id VARCHAR NOT NULL,
            PRIMARY KEY (rhea_id, pubmed_id)
        )
        """
    )
    _insert_rows_batched(
        connection,
        "reaction_publication",
        [
            (reaction.rhea_id, citation)
            for reaction in snapshot.reactions.values()
            for citation in sorted(reaction.citations)
        ],
    )


def _load_directions(
    connection: duckdb.DuckDBPyConnection, file_directions: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE reaction_quartet AS
        SELECT
            CAST(RHEA_ID_MASTER AS BIGINT) AS master_id,
            CAST(RHEA_ID_LR AS BIGINT) AS rhea_id_lr,
            CAST(RHEA_ID_RL AS BIGINT) AS rhea_id_rl,
            CAST(RHEA_ID_BI AS BIGINT) AS rhea_id_bi
        FROM read_csv(
            ?, delim = '\t', header = true, all_varchar = true,
            compression = ?
        )
        """,
        [str(file_directions), _duckdb_compression(file_directions)],
    )
    connection.execute("ALTER TABLE reaction_quartet ADD PRIMARY KEY (master_id)")


def _validate_rdf_directions(connection: duckdb.DuckDBPyConnection) -> None:
    mismatch_count = _fetch_scalar_int(
        connection,
        """
        WITH expected AS (
            SELECT master_id AS rhea_id, master_id, 'UN' AS direction
            FROM reaction_quartet
            UNION ALL
            SELECT rhea_id_lr, master_id, 'LR' FROM reaction_quartet
            UNION ALL
            SELECT rhea_id_rl, master_id, 'RL' FROM reaction_quartet
            UNION ALL
            SELECT rhea_id_bi, master_id, 'BI' FROM reaction_quartet
        )
        SELECT count(*)
        FROM expected
        LEFT JOIN reaction USING (rhea_id)
        WHERE reaction.rhea_id IS NULL
           OR reaction.master_id != expected.master_id
           OR reaction.direction != expected.direction
        """,
    )
    if mismatch_count:
        raise ValueError(
            "Rhea RDF direction semantics disagree with "
            f"rhea-directions.tsv ({mismatch_count} rows)"
        )


def _load_relationships(
    connection: duckdb.DuckDBPyConnection, file_relationships: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE reaction_relationship AS
        SELECT
            CAST(FROM_REACTION_ID AS BIGINT) AS from_reaction_id,
            CAST(TO_REACTION_ID AS BIGINT) AS to_reaction_id,
            TYPE AS relation_type
        FROM read_csv(
            ?, delim = '\t', header = true, all_varchar = true,
            compression = ?
        )
        """,
        [str(file_relationships), _duckdb_compression(file_relationships)],
    )


def _load_obsoletes(
    connection: duckdb.DuckDBPyConnection, file_obsoletes: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE obsolete_reaction AS
        SELECT CAST(RHEA_ID AS BIGINT) AS rhea_id
        FROM read_csv(
            ?, delim = '\t', header = true, all_varchar = true,
            compression = ?
        )
        """,
        [str(file_obsoletes), _duckdb_compression(file_obsoletes)],
    )
    if _table_exists(connection, "reaction"):
        connection.execute(
            """
            UPDATE reaction
            SET is_obsolete = true
            WHERE rhea_id IN (SELECT rhea_id FROM obsolete_reaction)
            """
        )


def _load_reaction_smiles(
    connection: duckdb.DuckDBPyConnection, file_smiles: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE reaction_smiles AS
        SELECT rhea_id, reaction_smiles
        FROM read_csv(
            ?,
            delim = '\t',
            header = false,
            compression = ?,
            columns = {
                'rhea_id': 'BIGINT',
                'reaction_smiles': 'VARCHAR'
            }
        )
        """,
        [str(file_smiles), _duckdb_compression(file_smiles)],
    )


def _load_sdf(connection: duckdb.DuckDBPyConnection, file_sdf: Path) -> None:
    connection.execute(
        """
        CREATE TABLE compound_structure (
            accession VARCHAR PRIMARY KEY,
            role VARCHAR,
            chebi_xref VARCHAR,
            generic_compound_accession VARCHAR,
            underlying_chebi_polymer_accession VARCHAR,
            formula VARCHAR,
            charge_text VARCHAR,
            charge_numeric INTEGER,
            name VARCHAR,
            molblock VARCHAR NOT NULL
        )
        """
    )
    _insert_rows_batched(
        connection,
        "compound_structure",
        [
            (
                structure.accession,
                structure.role,
                structure.chebi_xref,
                structure.generic_compound_accession,
                structure.underlying_chebi_polymer_accession,
                structure.formula,
                structure.charge,
                _optional_int(structure.charge),
                structure.name,
                structure.molblock,
            )
            for structure in iter_sdf_structures(file_sdf)
        ],
    )


def _load_chebi_names(
    connection: duckdb.DuckDBPyConnection, file_chebi_names: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE chebi_name AS
        SELECT
            'CHEBI:' || CAST(
                CAST(replace(upper(trim(chebi_id)), 'CHEBI:', '') AS BIGINT)
                AS VARCHAR
            ) AS chebi_id,
            trim(name) AS name
        FROM read_csv(
            ?,
            delim = '\t',
            header = false,
            compression = ?,
            columns = {'chebi_id': 'VARCHAR', 'name': 'VARCHAR'}
        )
        """,
        [str(file_chebi_names), _duckdb_compression(file_chebi_names)],
    )


def _load_chebi_mapping(
    connection: duckdb.DuckDBPyConnection, file_mapping: Path
) -> None:
    connection.execute(
        """
        CREATE TABLE chebi_ph7_3_mapping AS
        SELECT
            'CHEBI:' || CAST(CAST(CHEBI AS BIGINT) AS VARCHAR) AS chebi_id,
            'CHEBI:' || CAST(CAST(CHEBI_PH7_3 AS BIGINT) AS VARCHAR)
                AS chebi_ph7_3_id,
            ORIGIN AS origin
        FROM read_csv(
            ?, delim = '\t', header = true, all_varchar = true,
            compression = ?
        )
        """,
        [str(file_mapping), _duckdb_compression(file_mapping)],
    )


def _load_xrefs(connection: duckdb.DuckDBPyConnection, file_xrefs: Path) -> None:
    connection.execute(
        """
        CREATE TABLE reaction_xref AS
        SELECT
            CAST(RHEA_ID AS BIGINT) AS rhea_id,
            DIRECTION AS direction,
            CAST(MASTER_ID AS BIGINT) AS master_id,
            ID AS external_id,
            DB AS database_name
        FROM read_csv(
            ?, delim = '\t', header = true, all_varchar = true,
            compression = ?
        )
        """,
        [str(file_xrefs), _duckdb_compression(file_xrefs)],
    )


def _load_uniprot(
    connection: duckdb.DuckDBPyConnection,
    *,
    file_sprot: Path | None,
    file_trembl: Path | None,
) -> None:
    connection.execute(
        """
        CREATE TABLE reaction_uniprot (
            rhea_id BIGINT NOT NULL,
            direction VARCHAR NOT NULL,
            master_id BIGINT NOT NULL,
            uniprot_id VARCHAR NOT NULL,
            uniprot_section VARCHAR NOT NULL
        )
        """
    )
    for path, section in (
        (file_sprot, "Swiss-Prot"),
        (file_trembl, "TrEMBL"),
    ):
        if path is None:
            continue
        connection.execute(
            """
            INSERT INTO reaction_uniprot
            SELECT
                CAST(RHEA_ID AS BIGINT),
                DIRECTION,
                CAST(MASTER_ID AS BIGINT),
                ID,
                ?
            FROM read_csv(
                ?, delim = '\t', header = true, all_varchar = true,
                compression = ?
            )
            """,
            [section, str(path), _duckdb_compression(path)],
        )


def _validate_specialized_xref(
    connection: duckdb.DuckDBPyConnection,
    file_xref: Path,
    *,
    database_name: str,
) -> None:
    mismatch_count = _fetch_scalar_int(
        connection,
        """
        WITH specialized AS (
            SELECT
                CAST(RHEA_ID AS BIGINT) AS rhea_id,
                DIRECTION AS direction,
                CAST(MASTER_ID AS BIGINT) AS master_id,
                ID AS external_id
            FROM read_csv(
                ?, delim = '\t', header = true, all_varchar = true,
                compression = ?
            )
        ),
        aggregate_xref AS (
            SELECT rhea_id, direction, master_id, external_id
            FROM reaction_xref
            WHERE database_name = ?
        ),
        difference AS (
            (SELECT * FROM specialized EXCEPT SELECT * FROM aggregate_xref)
            UNION ALL
            (SELECT * FROM aggregate_xref EXCEPT SELECT * FROM specialized)
        )
        SELECT count(*) FROM difference
        """,
        [str(file_xref), _duckdb_compression(file_xref), database_name],
    )
    if mismatch_count:
        raise ValueError(
            f"Rhea {database_name} specialized mapping disagrees with "
            f"rhea2xrefs.tsv ({mismatch_count} rows)"
        )


def _create_views(connection: duckdb.DuckDBPyConnection) -> None:
    if _table_exists(connection, "reaction_xref"):
        connection.execute(
            """
            CREATE VIEW reaction_ec AS
            SELECT rhea_id, direction, master_id, external_id AS ec_number
            FROM reaction_xref
            WHERE database_name = 'EC'
            """
        )
        connection.execute(
            """
            CREATE VIEW reaction_go AS
            SELECT rhea_id, direction, master_id, external_id AS go_id
            FROM reaction_xref
            WHERE database_name = 'GO'
            """
        )
    if _table_exists(connection, "reaction") and _table_exists(
        connection, "reaction_participant"
    ):
        connection.execute(
            """
            CREATE VIEW reaction_participant_direction AS
            SELECT
                reaction.rhea_id,
                reaction.master_id,
                reaction.direction,
                participant.participant_id,
                participant.compound_id,
                participant.side,
                participant.coefficient_text,
                participant.coefficient_numeric,
                participant.location,
                CASE
                    WHEN reaction.direction = 'LR' AND participant.side = 'L'
                        THEN 'substrate'
                    WHEN reaction.direction = 'LR' AND participant.side = 'R'
                        THEN 'product'
                    WHEN reaction.direction = 'RL' AND participant.side = 'R'
                        THEN 'substrate'
                    WHEN reaction.direction = 'RL' AND participant.side = 'L'
                        THEN 'product'
                    ELSE NULL
                END AS directional_role
            FROM reaction
            JOIN reaction_participant AS participant
              ON participant.master_id = reaction.master_id
            """
        )


def _record_metadata(
    connection: duckdb.DuckDBPyConnection,
    *,
    sources: Mapping[str, Path],
    display_paths: Mapping[str, str],
    scope: str,
    release_number: int | None,
    release_date: str | None,
) -> None:
    try:
        bioextract_version = metadata.version("bioextract")
    except metadata.PackageNotFoundError:
        bioextract_version = "unknown"
    source_rows = [
        (
            name,
            display_paths[name],
            None,
            infer_media_type(path),
            None,
        )
        for name, path in sources.items()
    ]
    values = {
        "bioextract.metadata_schema_version": METADATA_SCHEMA_VERSION,
        "bioextract.resource_name": "rhea",
        "bioextract.resource_schema_version": SCHEMA_VERSION,
        "bioextract.source_schema_profile": SOURCE_SCHEMA_PROFILE,
        "bioextract.scope": scope,
        "bioextract.generated_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "bioextract.package_version": bioextract_version,
        "bioextract.release_version": (
            None if release_number is None else str(release_number)
        ),
        "bioextract.release_version_source": (
            None if release_number is None else "official_metadata"
        ),
        "bioextract.rhea_release_date": release_date,
        "bioextract.validation_status": "passed",
        "bioextract.validation_issue_count": "0",
    }
    connection.executemany(
        "INSERT INTO _bioextract.metadata VALUES (?, ?)",
        [(key, value) for key, value in values.items() if value is not None],
    )
    connection.executemany(
        "INSERT INTO _bioextract.source_file VALUES (?, ?, ?, ?, ?)",
        source_rows,
    )


def _record_table_inventory(
    connection: duckdb.DuckDBPyConnection,
) -> dict[str, int]:
    names = [
        row[0]
        for row in connection.execute(
            """
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'main' AND table_type = 'BASE TABLE'
            ORDER BY table_name
            """
        ).fetchall()
    ]
    counts = {
        name: _fetch_scalar_int(connection, f'SELECT count(*) FROM "{name}"')
        for name in names
    }
    connection.executemany(
        "INSERT INTO _bioextract.table_info VALUES (?, ?, ?)",
        [(name, "canonical", count) for name, count in counts.items()],
    )
    return counts


def _validate_staged_database(
    path: Path,
    *,
    row_counts: Mapping[str, int],
) -> None:
    connection = _connect_publication(path, read_only=True)
    try:
        connection.execute("PRAGMA database_size").fetchall()
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        validate_duckdb_metadata_v2(connection, metadata)
        required_metadata = {
            "bioextract.metadata_schema_version",
            "bioextract.resource_name",
            "bioextract.resource_schema_version",
            "bioextract.source_schema_profile",
            "bioextract.scope",
            "bioextract.generated_at",
            "bioextract.package_version",
        }
        missing = sorted(required_metadata - set(metadata))
        if missing:
            raise RuntimeError(
                f"Rhea publication is missing provenance keys: {missing}"
            )
        inventory = dict(
            connection.execute(
                "SELECT table_name, row_count FROM _bioextract.table_info"
            ).fetchall()
        )
        if inventory != dict(row_counts):
            raise RuntimeError("Rhea table inventory does not match row counts")
        for table_name, expected_count in row_counts.items():
            observed = _fetch_scalar_int(
                connection,
                f'SELECT count(*) FROM "{table_name}"',
            )
            if observed != expected_count:
                raise RuntimeError(
                    f"Rhea row-count validation failed for table: {table_name}"
                )
    finally:
        connection.close()


def _connect_publication(
    path: Path,
    *,
    read_only: bool = False,
) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(
        str(path),
        read_only=read_only,
        config={"threads": str(pl.thread_pool_size())},
    )


def _release_metadata(
    sources: Mapping[str, Path],
) -> tuple[int | None, str | None]:
    file_properties = sources.get("release_properties")
    if file_properties is None:
        return None, None
    return read_release_properties(file_properties)


def _table_exists(connection: duckdb.DuckDBPyConnection, table_name: str) -> bool:
    return bool(
        _fetch_scalar_int(
            connection,
            """
            SELECT count(*)
            FROM information_schema.tables
            WHERE table_schema = 'main' AND table_name = ?
            """,
            [table_name],
        )
    )


def _optional_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _chebi_curie(value: int | str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.upper().startswith("CHEBI:"):
        text = text.split(":", maxsplit=1)[1]
    try:
        return f"CHEBI:{int(text)}"
    except ValueError:
        return None


def _duckdb_compression(file_path: Path) -> str:
    return "gzip" if is_gzip_file(file_path) else "none"


def _fetch_scalar_int(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    parameters: list[str] | None = None,
) -> int:
    row = connection.execute(query, parameters or []).fetchone()
    if row is None:
        raise RuntimeError("Expected one scalar query result")
    return int(row[0])


def _insert_rows_batched(
    connection: duckdb.DuckDBPyConnection,
    table_name: str,
    rows: Iterable[tuple[object, ...]],
) -> None:
    rows_materialized = list(rows)
    if not rows_materialized:
        return
    columns = [
        str(row[0])
        for row in connection.execute(
            """
            SELECT column_name
            FROM information_schema.columns
            WHERE table_schema = 'main' AND table_name = ?
            ORDER BY ordinal_position
            """,
            [table_name],
        ).fetchall()
    ]
    frame = pl.DataFrame(
        rows_materialized,
        schema=columns,
        orient="row",
        infer_schema_length=None,
    )
    with tempfile.NamedTemporaryFile(
        prefix="bioextract-rhea-batch-",
        suffix=".parquet",
        delete=False,
    ) as handle:
        file_batch = Path(handle.name)
    try:
        frame.write_parquet(file_batch)
        connection.execute(
            f'INSERT INTO "{table_name}" SELECT * FROM read_parquet(?)',
            [str(file_batch)],
        )
    finally:
        file_batch.unlink(missing_ok=True)


def _remove_staging_files(staging: Path) -> None:
    staging.unlink(missing_ok=True)
    Path(f"{staging}.wal").unlink(missing_ok=True)
