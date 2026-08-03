from pathlib import Path

import duckdb
import pytest

from bioextract import GODatabase
from bioextract.errors import CapabilityError, IntegrityError


def _write_publication(tmp_path: Path) -> Path:
    source = tmp_path / "go.obo"
    source.write_text(
        """format-version: 1.2

[Term]
id: GO:0000001
name: root process
namespace: biological_process
def: "root" [GOC:test]

[Term]
id: GO:0000002
name: child process
namespace: biological_process
def: "child" [GOC:test]
is_a: GO:0000001 ! root process
alt_id: GO:1234567
""",
        encoding="utf-8",
    )
    publication = tmp_path / "go.duckdb"
    GODatabase.from_obo(source).write_duckdb(publication)
    return publication


def test_source_handle_does_not_expose_native_connection(tmp_path: Path) -> None:
    source = tmp_path / "go.obo"
    source.write_text("format-version: 1.2\n", encoding="utf-8")

    with pytest.raises(CapabilityError, match="from_duckdb"):
        GODatabase.from_obo(source).connect()


@pytest.mark.parametrize(
    ("metadata_key", "value", "message"),
    [
        ("bioextract.resource_name", "kegg", "GO publication"),
        ("bioextract.source_schema_profile", "forged-v1", "source schema profile"),
        ("bioextract.resource_schema_version", "forged-v1", "schema version"),
    ],
)
def test_from_duckdb_rejects_wrong_go_identity_profile_and_schema(
    tmp_path: Path,
    metadata_key: str,
    value: str,
    message: str,
) -> None:
    publication = _write_publication(tmp_path)
    with duckdb.connect(str(publication)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? WHERE key=?",
            [value, metadata_key],
        )

    with pytest.raises(IntegrityError, match=message):
        GODatabase.from_duckdb(publication)


def test_from_duckdb_rejects_inventory_role_and_physical_schema_drift(
    tmp_path: Path,
) -> None:
    publication = _write_publication(tmp_path)
    with duckdb.connect(str(publication)) as connection:
        connection.execute("CREATE VIEW unrecorded AS SELECT * FROM term")

    with pytest.raises(IntegrityError, match="table/view inventory"):
        GODatabase.from_duckdb(publication)

    publication.unlink()
    publication = _write_publication(tmp_path)
    with duckdb.connect(str(publication)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info SET table_role='forged' "
            "WHERE table_name='term'"
        )
    with pytest.raises(IntegrityError, match="capability inventory"):
        GODatabase.from_duckdb(publication)

    publication.unlink()
    publication = _write_publication(tmp_path)
    with duckdb.connect(str(publication)) as connection:
        connection.execute("ALTER TABLE term ALTER go_id TYPE INTEGER USING 1")
    with pytest.raises(IntegrityError, match="table schema"):
        GODatabase.from_duckdb(publication)


def test_reopened_handle_rejects_replaced_publication(tmp_path: Path) -> None:
    publication = _write_publication(tmp_path)
    reopened = GODatabase.from_duckdb(publication)
    reopened.build_tidy()
    replacement_dir = tmp_path / "replacement"
    replacement_dir.mkdir()
    replacement = _write_publication(replacement_dir)
    replacement.replace(publication)

    with pytest.raises(IntegrityError, match="replaced"):
        reopened.connect()
    with pytest.raises(IntegrityError, match="replaced"):
        reopened.select_terms(term_ids=["GO:0000001"])
    with pytest.raises(IntegrityError, match="replaced"):
        reopened.build_tidy()


def test_reopened_handle_translates_vanished_publication_to_integrity_error(
    tmp_path: Path,
) -> None:
    publication = _write_publication(tmp_path)
    reopened = GODatabase.from_duckdb(publication)
    publication.unlink()

    with pytest.raises(IntegrityError, match="replaced"):
        reopened.connect()
