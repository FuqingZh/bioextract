from __future__ import annotations

import gzip
import zipfile
from pathlib import Path

import duckdb
import pytest

from bioextract import ChEBIDatabase
from bioextract.errors import IntegrityError

COMPOUNDS = """id\tname\tstatus_id\tsource\tparent_id\tmerge_type\tchebi_accession\tdefinition\tascii_name\tstars\tmodified_on\trelease_date
1\twater\t1\tChEBI\t\t\tCHEBI:1\twater definition\twater\t3\t2026-01-01\t2026-01-01
"""
NAMES = """id\tcompound_id\tname\ttype\tstatus_id\tadapted\tlanguage_code\tascii_name
10\t1\taqua\tSYNONYM\t1\tF\ten\taqua
"""
DOMAIN_OBO = """format-version: 1.2

[Term]
id: CHEBI:1
name: water
def: "water definition" []
alt_id: CHEBI:100
subset: 3_STAR
synonym: "aqua" EXACT []
xref: kegg.compound:C00001
xref: HMDB:HMDB0002111
property_value: http://purl.obolibrary.org/obo/chebi/formula "H2O" xsd:string
property_value: http://purl.obolibrary.org/obo/chebi/inchi "InChI=1S/H2O/h1H2" xsd:string
property_value: http://purl.obolibrary.org/obo/chebi/inchikey "XLYOFNOQVPJJNP-UHFFFAOYSA-N" xsd:string
is_a: CHEBI:2 ! parent

[Term]
id: CHEBI:2
name: parent
subset: 2_STAR

[Term]
id: CHEBI:3
name: obsolete
subset: 3_STAR
is_obsolete: true

[Term]
id: CHEBI:4
name: low star
subset: 1_STAR

[Term]
id: CHEBI:5
name: orphan source
subset: 3_STAR
alt_id: CHEBI:500
is_a: CHEBI:999 ! missing target
"""


def _write_release(directory: Path, *, gzip_compounds: bool = False) -> None:
    directory.mkdir()
    if gzip_compounds:
        with gzip.open(
            directory / "compounds.tsv.gz",
            "wt",
            encoding="utf-8",
        ) as handle:
            handle.write(COMPOUNDS)
    else:
        (directory / "compounds.tsv").write_text(COMPOUNDS, encoding="utf-8")
    (directory / "names.tsv").write_text(NAMES, encoding="utf-8")


def _write_domain_publication(tmp_path: Path) -> Path:
    file_obo = tmp_path / "chebi.obo"
    file_obo.write_text(DOMAIN_OBO, encoding="utf-8")
    file_sdf = tmp_path / "chebi.sdf"
    file_sdf.write_text(
        """water
  bioextract

  0  0  0  0  0  0            999 V2000
M  END
> <ChEBI ID>
CHEBI:1

$$$$
orphan
  bioextract

  0  0  0  0  0  0            999 V2000
M  END
> <ChEBI ID>
CHEBI:999

$$$$
""",
        encoding="utf-8",
    )
    path = tmp_path / "chebi.duckdb"
    ChEBIDatabase.from_obo(file_obo, sdf=file_sdf).write_duckdb(path)
    return path


def test_release_tables_keep_official_headers_and_use_internal_metadata(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "release"
    _write_release(directory, gzip_compounds=True)

    path = tmp_path / "chebi.duckdb"
    result = ChEBIDatabase.from_release(directory).write_duckdb(path)

    assert result.tables == ("compound", "compound_name")
    with duckdb.connect(str(path), read_only=True) as connection:
        assert [
            row[1]
            for row in connection.execute("PRAGMA table_info('compound')").fetchall()
        ] == [
            "id",
            "name",
            "status_id",
            "source",
            "parent_id",
            "merge_type",
            "chebi_accession",
            "definition",
            "ascii_name",
            "stars",
            "modified_on",
            "release_date",
        ]
        assert connection.execute(
            "SELECT table_role, row_count FROM _bioextract.table_info "
            "WHERE table_name = 'compound'"
        ).fetchone() == ("entity", 1)
        assert connection.execute(
            "SELECT count(*) FROM _bioextract.column_mapping"
        ).fetchone() == (0,)


def test_unsafe_release_archive_member_is_rejected(tmp_path: Path) -> None:
    file_archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(file_archive, "w") as archive:
        archive.writestr("../compounds.tsv", COMPOUNDS)
    with pytest.raises(ValueError, match="Unsafe archive member"):
        ChEBIDatabase.from_release(file_archive)


def test_validation_issue_metadata_and_canonical_fail_fast(tmp_path: Path) -> None:
    path = _write_domain_publication(tmp_path)
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key='bioextract.validation_status'"
        ).fetchone() == ("passed_with_warnings",)
        assert connection.execute(
            "SELECT issue_code, relation_name FROM _bioextract.validation_issue "
            "ORDER BY issue_id"
        ).fetchall() == [
            ("foreign_key_violation", "compound_relation"),
            ("foreign_key_violation", "compound_structure"),
        ]

    invalid = tmp_path / "invalid.obo"
    invalid.write_text(
        """[Term]
id: CHEBI:1
name: first

[Term]
id: CHEBI:1
name: duplicate
""",
        encoding="utf-8",
    )
    preserved = tmp_path / "preserved.duckdb"
    preserved.write_bytes(b"old")
    with pytest.raises(IntegrityError, match="duplicate"):
        ChEBIDatabase.from_obo(invalid).write_duckdb(
            preserved,
            if_exists="replace",
        )
    assert preserved.read_bytes() == b"old"
    assert not list(tmp_path.glob(".preserved.duckdb.*"))


def test_chebi_rejects_legacy_and_unknown_metadata_contracts(
    tmp_path: Path,
) -> None:
    current = _write_domain_publication(tmp_path)
    legacy = tmp_path / "legacy.duckdb"
    legacy.write_bytes(current.read_bytes())
    with duckdb.connect(str(legacy)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='1' "
            "WHERE key='bioextract.metadata_schema_version'"
        )
        connection.execute(
            "UPDATE _bioextract.metadata SET key='bioextract.schema_version' "
            "WHERE key='bioextract.resource_schema_version'"
        )
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="five _bioextract relations"):
        ChEBIDatabase.from_duckdb(legacy)

    unknown = tmp_path / "unknown.duckdb"
    unknown.write_bytes(current.read_bytes())
    with duckdb.connect(str(unknown)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='3' "
            "WHERE key='bioextract.metadata_schema_version'"
        )
    with pytest.raises(ValueError, match="metadata schema version"):
        ChEBIDatabase.from_duckdb(unknown)

    missing_v1_table = tmp_path / "missing-v1-table.duckdb"
    missing_v1_table.write_bytes(current.read_bytes())
    with duckdb.connect(str(missing_v1_table)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="five _bioextract relations"):
        ChEBIDatabase.from_duckdb(missing_v1_table)

    for index, profile in enumerate(("", "unknown-profile")):
        unsupported_profile = tmp_path / f"unsupported-profile-{index}.duckdb"
        unsupported_profile.write_bytes(current.read_bytes())
        with duckdb.connect(str(unsupported_profile)) as connection:
            connection.execute(
                "UPDATE _bioextract.metadata SET value=? "
                "WHERE key='bioextract.source_schema_profile'",
                [profile],
            )
        with pytest.raises(ValueError, match="source schema profile"):
            ChEBIDatabase.from_duckdb(unsupported_profile)
