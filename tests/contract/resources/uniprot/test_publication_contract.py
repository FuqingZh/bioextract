from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.uniprot import UniProtDatabase


def _write_dat(path: Path) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            """ID   TEST_HUMAN Reviewed; 3 AA.
AC   P12345;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
     ACD
//
"""
        )
    return path


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for identifier, sequence in records.items():
            handle.write(f">sp|{identifier}|fixture\n{sequence}\n")
    return path


def _write_idmapping(path: Path) -> Path:
    row = [
        "P12345",
        "TEST_HUMAN",
        "1234",
        "NP_000001.1",
        "",
        "",
        "GO:0003677",
        "",
        "",
        "",
        "",
        "",
        "9606",
        "",
        "",
        "",
        "",
        "",
        "ENSG1",
        "ENST1",
        "ENSP1",
        "",
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(row) + "\n")
    return path


def test_idmapping_schema_and_path_validation(tmp_path: Path) -> None:
    bad = tmp_path / "bad.parquet"
    pl.DataFrame({"UniProtId": ["P12345"]}).write_parquet(bad)
    with pytest.raises(ValueError, match="schema mismatch"):
        UniProtDatabase.from_idmapping(bad).scan_mapping()
    with pytest.raises(FileNotFoundError):
        UniProtDatabase.from_idmapping(tmp_path / "missing.tab.gz")
    unsupported = tmp_path / "mapping.txt"
    unsupported.write_text("", encoding="utf-8")
    with pytest.raises(pl.exceptions.NoDataError, match="empty CSV"):
        UniProtDatabase.from_idmapping(unsupported).scan_mapping().collect_schema()
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no parquet"):
        UniProtDatabase.from_idmapping(empty)


@pytest.mark.parametrize("container", ["parquet", "hive"])
def test_idmapping_parquet_schema_requires_all_string_columns(
    tmp_path: Path, container: str
) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    frame = (
        UniProtDatabase.from_idmapping(raw)
        .scan_mapping()
        .collect()
        .with_columns(pl.col("GeneId").cast(pl.Int64))
    )
    if container == "parquet":
        source = tmp_path / "wrong-type.parquet"
        frame.write_parquet(source)
    else:
        source = tmp_path / "wrong-type-hive"
        partition = source / "TaxId=9606"
        partition.mkdir(parents=True)
        frame.drop("TaxId").write_parquet(partition / "part.parquet")

    with pytest.raises(ValueError, match="schema mismatch"):
        UniProtDatabase.from_idmapping(source).scan_mapping()


def test_idmapping_publication_is_atomic(tmp_path: Path) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    destination = tmp_path / "mapping.parquet"
    destination.write_bytes(b"existing")
    database = UniProtDatabase.from_idmapping(raw)
    with pytest.raises(FileExistsError):
        database.write_parquet(destination, taxon_ids=["9606"])
    assert destination.read_bytes() == b"existing"


def test_roles_are_validated_by_content_not_basename(tmp_path: Path) -> None:
    entries = _write_dat(tmp_path / "anything.gz")
    canonical = _write_fasta(tmp_path / "other.dat", {"P12345": "WRONG"})
    with pytest.raises(ValueError, match="Canonical FASTA"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, canonical_sequences=canonical
        ).write_duckdb(tmp_path / "bad.duckdb")

    with pytest.raises(ValueError, match="no UniProtKB records"):
        UniProtDatabase.from_knowledgebase(
            entries=_write_fasta(tmp_path / "not-dat", {"P12345": "ACD"})
        ).write_duckdb(tmp_path / "wrong-role.duckdb")


def test_release_version_is_never_inferred_from_path(tmp_path: Path) -> None:
    entries = _write_dat(tmp_path / "2026_01-uniprot.dat.gz")
    path = tmp_path / "unknown-version.duckdb"
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(path)
    assert UniProtDatabase.from_duckdb(path).release_version is None


def test_missing_required_source_field_always_fails(tmp_path: Path) -> None:
    path = tmp_path / "missing-ox.dat"
    path.write_text(
        """ID   TEST_HUMAN Reviewed; 3 AA.
AC   P12345;
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
     ACD
//
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="OX"):
        UniProtDatabase.from_knowledgebase(entries=path).write_duckdb(
            tmp_path / "invalid.duckdb"
        )


def test_existing_destination_fails_before_entries_are_parsed(tmp_path: Path) -> None:
    entries = tmp_path / "invalid.dat"
    entries.write_text("not a UniProt record\n", encoding="utf-8")
    destination = tmp_path / "existing.duckdb"
    destination.write_bytes(b"existing")

    with pytest.raises(FileExistsError):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(destination)
    assert destination.read_bytes() == b"existing"


def test_primary_accessions_must_be_unique_across_records(tmp_path: Path) -> None:
    entries = tmp_path / "duplicate-primary-accession.dat"
    entries.write_text(
        """ID   FIRST_HUMAN Reviewed; 3 AA.
AC   P11111;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
     ACD
//
ID   SECOND_HUMAN Reviewed; 3 AA.
AC   P11111;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  365 MW;  69CB1DB000000000 CRC64;
     AEF
//
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate primary UniProt accession: P11111"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "duplicate.duckdb"
        )


def test_duplicate_accession_within_record_is_rejected(tmp_path: Path) -> None:
    entries = tmp_path / "duplicate-record-accession.dat"
    entries.write_text(
        """ID   FIRST_HUMAN Reviewed; 3 AA.
AC   P11111; Q99999; Q99999;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
     ACD
//
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate accession in UniProtKB record 1"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "duplicate.duckdb"
        )


def test_metadata_v3_requires_validation_issue_table(tmp_path: Path) -> None:
    path = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="validation_issue"):
        UniProtDatabase.from_duckdb(path)


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        (
            "ALTER TABLE protein ALTER sequence_length TYPE VARCHAR",
            "physical schema mismatch",
        ),
        (
            "ALTER TABLE protein RENAME sequence_length TO reported_length",
            "physical schema mismatch",
        ),
        ("CREATE TABLE unexpected(value INTEGER)", "relation inventory mismatch"),
        (
            "INSERT INTO _bioextract.table_info VALUES ('unexpected', 'canonical', 0)",
            "table_info inventory mismatch",
        ),
        (
            "UPDATE _bioextract.metadata SET value='unsupported-profile' "
            "WHERE key='bioextract.source_schema_profile'",
            "source schema profile",
        ),
    ],
)
def test_from_duckdb_rejects_physical_contract_corruption(
    tmp_path: Path, corruption: str, message: str
) -> None:
    path = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(corruption)

    with pytest.raises(ValueError, match=message):
        UniProtDatabase.from_duckdb(path)
