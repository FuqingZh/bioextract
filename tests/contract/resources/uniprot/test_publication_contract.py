from __future__ import annotations

import gc
import gzip
import json
import os
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.errors import CapabilityError, IntegrityError
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


def test_idmapping_relative_source_is_pinned_for_lazy_reads(
    tmp_path: Path,
) -> None:
    source = _write_idmapping(tmp_path / "mapping.tab.gz")
    previous = Path.cwd()
    try:
        os.chdir(tmp_path)
        database = UniProtDatabase.from_idmapping(source.name)
        os.chdir(previous)
        assert database.read_mapping(taxon_ids=["9606"]).height == 1
    finally:
        os.chdir(previous)


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
    destination = tmp_path / "mapping.duckdb"
    destination.write_bytes(b"existing")
    database = UniProtDatabase.from_idmapping(raw)
    with pytest.raises(FileExistsError):
        database.write_duckdb(destination, taxon_ids=["9606"])
    assert destination.read_bytes() == b"existing"


def test_idmapping_duckdb_contract_and_reopen_parity(tmp_path: Path) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    source = UniProtDatabase.from_idmapping(raw, release_version="2026_01")
    path = tmp_path / "mapping.duckdb"
    result = source.write_duckdb(path, allow_all_taxa=True)

    assert result.tables == ("mapping",)
    assert not path.with_suffix(".json").exists()
    reopened = UniProtDatabase.from_duckdb(path)
    assert reopened.release_version == "2026_01"
    assert reopened.read_mapping(taxon_ids=["9606"]).equals(
        source.read_mapping(taxon_ids=["9606"])
    )
    with pytest.raises(ValueError, match="allow_all_taxa"):
        reopened.read_mapping()
    assert reopened.read_mapping(allow_all_taxa=True).height == 1
    first = reopened.connect()
    second = reopened.connect()
    try:
        assert first is not second
        assert first.execute("SELECT count(*) FROM mapping").fetchone() == (1,)
        metadata = dict(
            first.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "1"
        assert metadata["bioextract.source_schema_profile"] == (
            "uniprot-idmapping-selected-22-column-v1"
        )
        assert metadata["bioextract.capability.mapping"] == "true"
        assert json.loads(metadata["bioextract.scope"]) == {"all_taxa": True}
        assert {
            row[0]
            for row in first.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema='_bioextract'"
            ).fetchall()
        } == {
            "metadata",
            "source_file",
            "table_info",
            "column_mapping",
            "validation_issue",
        }
        with pytest.raises(duckdb.Error):
            first.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        first.close()
        second.close()


def test_from_duckdb_discriminates_uniprot_profiles(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='uniprotkb-flat-file-v1' "
            "WHERE key='bioextract.source_schema_profile'"
        )
    with pytest.raises(IntegrityError, match="source schema profile"):
        UniProtDatabase.from_duckdb(path)


def test_profile_capabilities_and_selection_honor_pinned_identity(
    tmp_path: Path,
) -> None:
    mapping_path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(mapping_path, allow_all_taxa=True)
    mapping = UniProtDatabase.from_duckdb(mapping_path)
    with pytest.raises(CapabilityError):
        mapping.select_ids(["P12345"], namespace="uniprot")

    knowledgebase_path = tmp_path / "knowledgebase.duckdb"
    replacement = tmp_path / "replacement.duckdb"
    source = UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    )
    source.write_duckdb(knowledgebase_path)
    source.write_duckdb(replacement)
    knowledgebase = UniProtDatabase.from_duckdb(knowledgebase_path)
    os.replace(replacement, knowledgebase_path)

    with pytest.raises(IntegrityError, match="replaced"):
        knowledgebase.select_ids(["P12345"], namespace="uniprot").extract_proteins()

    vanished = UniProtDatabase.from_duckdb(mapping_path)
    mapping_path.unlink()
    with pytest.raises(IntegrityError, match="unavailable"):
        vanished.connect()


def test_idmapping_reopen_requires_exact_source_role(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.source_file SET logical_name='fabricated_mapping'"
        )
        connection.execute(
            "UPDATE _bioextract.metadata SET value=replace(value, ?, ?) "
            "WHERE key='bioextract.sources'",
            ["idmapping_selected", "fabricated_mapping"],
        )
    with pytest.raises(IntegrityError, match="source inventory"):
        UniProtDatabase.from_duckdb(path)


def test_idmapping_reopen_rejects_fabricated_column_provenance(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "INSERT INTO _bioextract.column_mapping VALUES "
            "('mapping', 'UniProtId', 'fabricated', 'forged')"
        )
    with pytest.raises(IntegrityError, match="column provenance"):
        UniProtDatabase.from_duckdb(path)


def test_idmapping_publication_records_and_validates_selected_taxa(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, taxon_ids=[9606, "10090"])
    with duckdb.connect(str(path)) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert json.loads(metadata["bioextract.scope"]) == {
            "taxon_ids": ["9606", "10090"]
        }
        connection.execute(
            "UPDATE _bioextract.metadata SET value='{}' WHERE key='bioextract.scope'"
        )
    with pytest.raises(IntegrityError, match="taxon scope"):
        UniProtDatabase.from_duckdb(path)


def test_idmapping_reopen_rejects_non_normalized_taxon_scope(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? WHERE key='bioextract.scope'",
            [json.dumps({"taxon_ids": [" "]})],
        )
    with pytest.raises(IntegrityError, match="taxon scope"):
        UniProtDatabase.from_duckdb(path)


def test_idmapping_reopen_translates_malformed_metadata_schema(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute("ALTER TABLE _bioextract.metadata RENAME key TO malformed")
    with pytest.raises(IntegrityError):
        UniProtDatabase.from_duckdb(path)


@pytest.mark.parametrize("scope", [{"all_taxa": 1}, {"all_taxa": 1.0}])
def test_idmapping_reopen_requires_boolean_all_taxa(
    tmp_path: Path, scope: dict[str, object]
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? WHERE key='bioextract.scope'",
            [json.dumps(scope)],
        )
    with pytest.raises(IntegrityError, match="taxon scope"):
        UniProtDatabase.from_duckdb(path)


@pytest.mark.parametrize("sources", [{}, [{}]])
def test_idmapping_reopen_translates_malformed_embedded_sources(
    tmp_path: Path, sources: object
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? WHERE key='bioextract.sources'",
            [json.dumps(sources)],
        )
    with pytest.raises(IntegrityError):
        UniProtDatabase.from_duckdb(path)


def test_reopened_mapping_lazy_frame_owns_and_releases_connection(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    ).write_duckdb(path, allow_all_taxa=True)
    reopened = UniProtDatabase.from_duckdb(path)

    frame = reopened.scan_mapping(taxon_ids=["9606"])
    assert "PYTHON SCAN" in frame.explain()
    assert frame.collect().height == 1
    assert frame.collect().height == 1
    del frame
    gc.collect()

    with duckdb.connect(str(path)) as connection:
        connection.execute("CHECKPOINT")


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


def test_metadata_v1_requires_validation_issue_table(tmp_path: Path) -> None:
    path = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(IntegrityError, match="validation_issue"):
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
        (
            "UPDATE _bioextract.table_info SET table_role='derived' "
            "WHERE table_name='protein'",
            "role inventory",
        ),
        (
            "INSERT INTO _bioextract.metadata VALUES "
            "('bioextract.capability.unexpected', 'true')",
            "capability inventory",
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

    with pytest.raises(IntegrityError, match=message):
        UniProtDatabase.from_duckdb(path)
