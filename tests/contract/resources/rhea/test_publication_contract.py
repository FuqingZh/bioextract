from __future__ import annotations

import inspect
import os
import zipfile
from pathlib import Path
from typing import Any

import duckdb
import polars as pl
import pytest

from bioextract import RheaDatabase
from bioextract.errors import CapabilityError
from bioextract.rhea import _duckdb

RDF_FIXTURE = """\
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
  xmlns:rh="http://rdf.rhea-db.org/">
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000">
    <rh:accession>RHEA:10000</rh:accession>
    <rh:equation>A = B</rh:equation>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
    <rh:isChemicallyBalanced>true</rh:isChemicallyBalanced>
    <rh:isTransport>false</rh:isTransport>
    <rh:side rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:side rdf:resource="http://rdf.rhea-db.org/10000_R"/>
    <rh:citation rdf:resource="http://rdf.rhea-db.org/12345"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10001">
    <rh:accession>RHEA:10001</rh:accession>
    <rh:substrates rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:products rdf:resource="http://rdf.rhea-db.org/10000_R"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10002">
    <rh:accession>RHEA:10002</rh:accession>
    <rh:substrates rdf:resource="http://rdf.rhea-db.org/10000_R"/>
    <rh:products rdf:resource="http://rdf.rhea-db.org/10000_L"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10003">
    <rh:accession>RHEA:10003</rh:accession>
    <rh:substratesOrProducts rdf:resource="http://rdf.rhea-db.org/10000_L"/>
    <rh:substratesOrProducts rdf:resource="http://rdf.rhea-db.org/10000_R"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10004">
    <rh:accession>RHEA:10004</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Obsolete"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_L">
    <rh:curatedOrder>1</rh:curatedOrder>
    <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_A"/>
    <rh:contains1 rdf:resource="http://rdf.rhea-db.org/Participant_A"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_R">
    <rh:curatedOrder>2</rh:curatedOrder>
    <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_B"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/contains1">
    <rh:coefficient>2</rh:coefficient>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_A">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_A"/>
    <rh:location rdf:resource="http://rdf.rhea-db.org/Cytoplasm"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_B">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_B"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_A">
    <rh:id>1</rh:id>
    <rh:accession>CHEBI:1</rh:accession>
    <rh:name>A</rh:name>
    <rh:formula>C</rh:formula>
    <rh:charge>-1</rh:charge>
    <rh:chebi rdf:resource="http://purl.obolibrary.org/obo/CHEBI_1"/>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/SmallMolecule"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_B">
    <rh:id>2</rh:id>
    <rh:accession>GENERIC:2</rh:accession>
    <rh:name>B</rh:name>
    <rh:charge>(-1)n</rh:charge>
    <rh:reactivePart rdf:resource="http://rdf.rhea-db.org/Compound_Part"/>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/Polymer"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_Part">
    <rh:position>C1</rh:position>
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/ReactivePart"/>
  </rdf:Description>
</rdf:RDF>
"""

SDF_FIXTURE = """\
CHEBI:1
  bioextract

  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <ACCESSION>
CHEBI:1

> <ROLE>
small molecule

> <CHEBI_XREF>
CHEBI:1

> <Formula>
C

> <Charge>
(-1)n

> <Rhea_ascii_name>
A

$$$$
"""


def _write(path: Path, value: str) -> Path:
    path.write_text(value, encoding="utf-8")
    return path


def _write_release(root: Path) -> Path:
    raw = root / "raw"
    raw.mkdir(parents=True)
    _write(raw / "LICENSE.txt", "Creative Commons\n")
    _write(raw / "chebiId_name.tsv", "CHEBI:1\t A\n")
    _write(
        raw / "chebi_pH7_3_mapping.tsv",
        "CHEBI\tCHEBI_PH7_3\tORIGIN\n1\t1\tcomputation\n",
    )
    _write(
        raw / "rhea-directions.tsv",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "10000\t10001\t10002\t10003\n",
    )
    _write(raw / "rhea-obsoletes.tsv", "RHEA_ID\n")
    _write(raw / "rhea-reaction-smiles.tsv", "10000\tA>>B\n")
    _write(
        raw / "rhea-relationships.tsv",
        "FROM_REACTION_ID\tTYPE\tTO_REACTION_ID\n10000\tis_a\t10004\n",
    )
    _write(
        raw / "rhea-release.properties",
        "rhea.release.number=141\nrhea.release.date=2026-06-10\n",
    )
    _write(raw / "rhea.rdf", RDF_FIXTURE)
    _write(raw / "rhea.sdf", SDF_FIXTURE)
    header = "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n"
    _write(raw / "rhea2ec.tsv", header + "10000\tUN\t10000\t1.1.1.1\n")
    _write(raw / "rhea2go.tsv", header + "10000\tUN\t10000\tGO:0000001\n")
    _write(
        raw / "rhea2uniprot_sprot.tsv",
        header + "10000\tUN\t10000\tP00001\n",
    )
    _write(
        raw / "rhea2uniprot_trembl.tsv",
        header + "10001\tLR\t10000\tA00001\n",
    )
    _write(
        raw / "rhea2xrefs.tsv",
        "RHEA_ID\tDIRECTION\tMASTER_ID\tID\tDB\n"
        "10000\tUN\t10000\t1.1.1.1\tEC\n"
        "10000\tUN\t10000\tGO:0000001\tGO\n"
        "10000\tUN\t10000\tRXN-1\tECOCYC\n"
        "10000\tUN\t10000\tR00001\tKEGG_REACTION\n"
        "10000\tUN\t10000\tM0001\tMACIE\n"
        "10000\tUN\t10000\tMETA:RXN-1\tMETACYC\n"
        "10000\tUN\t10000\tR-HSA-1\tREACTOME\n",
    )
    return root


def test_from_files_signature_and_removed_constructors() -> None:
    parameters = inspect.signature(RheaDatabase.from_files).parameters
    assert tuple(parameters) == (
        "source",
        "rdf",
        "directions",
        "relationships",
        "obsolete_reactions",
        "reaction_smiles",
        "sdf",
        "chebi_names",
        "chebi_ph7_3_mapping",
        "xrefs",
        "uniprot_sprot",
        "uniprot_trembl",
    )
    assert parameters["source"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert parameters["source"].default is None
    assert all(
        parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is None
        for parameter in tuple(parameters.values())[1:]
    )
    assert not hasattr(RheaDatabase, "from_reaction_files")
    assert not hasattr(RheaDatabase, "from_compound_files")
    assert not hasattr(RheaDatabase, "from_cross_reference_files")
    assert not hasattr(RheaDatabase, "from_release")


def test_from_files_validation_and_scope_matrix(tmp_path: Path) -> None:
    files = {
        name: _write(tmp_path / f"{name}.txt", name)
        for name in (
            "rdf",
            "directions",
            "relationships",
            "obsolete_reactions",
            "reaction_smiles",
            "sdf",
            "chebi_names",
            "chebi_ph7_3_mapping",
            "xrefs",
            "uniprot_sprot",
            "uniprot_trembl",
        )
    }
    with pytest.raises(ValueError, match="At least one"):
        RheaDatabase.from_files()

    invalid_reaction_profiles = [
        {"rdf": files["rdf"]},
        {"directions": files["directions"]},
    ]
    for role in ("relationships", "obsolete_reactions", "reaction_smiles"):
        invalid_reaction_profiles.extend(
            (
                {role: files[role]},
                {"rdf": files["rdf"], role: files[role]},
                {"directions": files["directions"], role: files[role]},
            )
        )
    for values in invalid_reaction_profiles:
        with pytest.raises(ValueError, match="require both rdf and directions"):
            RheaDatabase.from_files(**values)

    reactions = {"rdf": files["rdf"], "directions": files["directions"]}
    assert RheaDatabase.from_files(**reactions).snapshot.scope == "reactions"
    for role in ("sdf", "chebi_names", "chebi_ph7_3_mapping"):
        assert RheaDatabase.from_files(**{role: files[role]}).snapshot.scope == (
            "compounds"
        )
    for role in ("xrefs", "uniprot_sprot", "uniprot_trembl"):
        assert RheaDatabase.from_files(**{role: files[role]}).snapshot.scope == (
            "cross_references"
        )

    compounds = {"chebi_names": files["chebi_names"]}
    cross_references = {"xrefs": files["xrefs"]}
    for values in (
        reactions | compounds,
        reactions | cross_references,
        compounds | cross_references,
        reactions | compounds | cross_references,
    ):
        assert RheaDatabase.from_files(**values).snapshot.scope == "partial"


def test_from_files_rejects_duplicate_physical_files(tmp_path: Path) -> None:
    source = _write(tmp_path / "source.tsv", "value\n")
    symlink = tmp_path / "alias.tsv"
    symlink.symlink_to(source)
    hard_link = tmp_path / "hard-link.tsv"
    os.link(source, hard_link)

    for alias in (tmp_path / "." / "source.tsv", symlink, hard_link):
        with pytest.raises(
            ValueError,
            match="roles 'chebi_names' and 'xrefs'.*same physical file",
        ):
            RheaDatabase.from_files(chebi_names=source, xrefs=alias)


def test_all_role_provenance_and_mixed_reaction_xref_capability(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    raw = release / "raw"
    path = tmp_path / "all-roles.duckdb"

    report = RheaDatabase.from_files(
        rdf=raw / "rhea.rdf",
        directions=raw / "rhea-directions.tsv",
        relationships=raw / "rhea-relationships.tsv",
        obsolete_reactions=raw / "rhea-obsoletes.tsv",
        reaction_smiles=raw / "rhea-reaction-smiles.tsv",
        sdf=raw / "rhea.sdf",
        chebi_names=raw / "chebiId_name.tsv",
        chebi_ph7_3_mapping=raw / "chebi_pH7_3_mapping.tsv",
        xrefs=raw / "rhea2xrefs.tsv",
        uniprot_sprot=raw / "rhea2uniprot_sprot.tsv",
        uniprot_trembl=raw / "rhea2uniprot_trembl.tsv",
    ).write_duckdb(path)

    assert report.scope == "partial"
    assert set(report.source_files) == {
        "rdf",
        "directions",
        "relationships",
        "obsoletes",
        "reaction_smiles",
        "sdf",
        "chebi_names",
        "chebi_ph7_3_mapping",
        "xrefs",
        "uniprot_sprot",
        "uniprot_trembl",
    }
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        logical_sources = {
            row[0]
            for row in connection.execute(
                "SELECT logical_name FROM _bioextract.source_file"
            ).fetchall()
        }
        assert metadata["bioextract.scope"] == "partial"
        assert metadata["bioextract.resource_schema_version"] == "rhea-duckdb-v1"
        assert "bioextract.release_version" not in metadata
        assert logical_sources == set(report.source_files)

    database = RheaDatabase.from_duckdb(path)
    assert database.snapshot.scope == "publication"
    matches = database.select_reactions(["1.1.1.1"], namespace="ec").matches().collect()
    assert matches["rhea_id"].to_list() == [10000]


def test_xref_only_publication_reports_reaction_capability_failure(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "xrefs-only.duckdb"
    RheaDatabase.from_files(xrefs=release / "raw" / "rhea2xrefs.tsv").write_duckdb(path)

    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    assert metadata["bioextract.scope"] == "cross_references"
    assert metadata["bioextract.resource_schema_version"] == "rhea-duckdb-v1"

    database = RheaDatabase.from_duckdb(path)
    assert database.snapshot.scope == "publication"
    with pytest.raises(CapabilityError, match="missing relations"):
        database.select_reactions(["1.1.1.1"], namespace="ec")


def test_rhea_publication_connections_share_the_polars_thread_budget(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    release = _write_release(tmp_path / "release")
    original_connect = duckdb.connect
    calls: list[tuple[bool, dict[str, Any] | None]] = []

    def connect(
        database: str | Path = ":memory:",
        read_only: bool = False,
        config: dict[str, Any] | None = None,
    ) -> duckdb.DuckDBPyConnection:
        calls.append((read_only, config))
        return original_connect(database=database, read_only=read_only, config=config)

    monkeypatch.setattr(_duckdb.duckdb, "connect", connect)

    RheaDatabase.from_files(xrefs=release / "raw" / "rhea2xrefs.tsv").write_duckdb(
        tmp_path / "rhea.duckdb"
    )

    expected_config = {"threads": str(pl.thread_pool_size())}
    assert calls == [(False, expected_config), (True, expected_config)]


def test_reaction_xref_publication_persists_partial_construction_scope(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    raw = release / "raw"
    path = tmp_path / "reaction-xrefs.duckdb"
    report = RheaDatabase.from_files(
        rdf=raw / "rhea.rdf",
        directions=raw / "rhea-directions.tsv",
        xrefs=raw / "rhea2xrefs.tsv",
    ).write_duckdb(path)

    assert report.scope == "partial"
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    assert metadata["bioextract.scope"] == "partial"
    assert metadata["bioextract.resource_schema_version"] == "rhea-duckdb-v1"

    database = RheaDatabase.from_duckdb(path)
    assert database.snapshot.scope == "publication"
    matches = database.select_reactions(["1.1.1.1"], namespace="ec").matches().collect()
    assert matches["rhea_id"].to_list() == [10000]


def test_release_is_strict_but_partial_constructor_is_not(
    tmp_path: Path,
) -> None:
    _write(tmp_path / "chebiId_name.tsv", "1\t A\n")

    with pytest.raises(ValueError, match="Incomplete Rhea release"):
        RheaDatabase.from_files(tmp_path)

    db = RheaDatabase.from_files(chebi_names=tmp_path / "chebiId_name.tsv")
    assert db.snapshot.scope == "compounds"


def test_source_backed_from_files_replaces_roles_and_records_final_provenance(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    overlay = _write(
        tmp_path / "overlay-directions.tsv",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "10000\t10001\t10002\t10003\n",
    )
    path = tmp_path / "source-overlay.duckdb"

    report = RheaDatabase.from_files(release, directions=overlay).write_duckdb(path)

    assert report.scope == "release"
    assert report.release_number == 141
    assert report.source_files["directions"] == str(overlay)
    assert str(release) in report.source_files["rdf"]
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    assert metadata["bioextract.scope"] == "release"


def test_source_backed_from_files_allows_only_overridden_missing_roles(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    missing = release / "raw" / "rhea.rdf"
    missing.unlink()
    overlay = _write(tmp_path / "overlay.rdf", RDF_FIXTURE)

    report = RheaDatabase.from_files(release, rdf=overlay).write_duckdb(
        tmp_path / "overridden-missing.duckdb"
    )
    assert report.scope == "release"
    assert report.source_files["rdf"] == str(overlay)

    missing_directions = _write_release(tmp_path / "missing-directions")
    (missing_directions / "raw" / "rhea-directions.tsv").unlink()
    with pytest.raises(ValueError, match="Incomplete Rhea release"):
        RheaDatabase.from_files(missing_directions)


def test_archive_source_overlay_replaces_missing_role_and_keeps_archive_provenance(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    archive_path = tmp_path / "release.zip"
    with zipfile.ZipFile(archive_path, mode="w") as archive:
        for path in release.rglob("*"):
            if path.is_file() and path.name != "rhea.rdf":
                archive.write(path, path.relative_to(tmp_path))
    overlay = _write(tmp_path / "overlay.rdf", RDF_FIXTURE)

    report = RheaDatabase.from_files(archive_path, rdf=overlay).write_duckdb(
        tmp_path / "archive-overlay.duckdb"
    )

    assert report.scope == "release"
    assert report.source_files["rdf"] == str(overlay)
    assert str(archive_path) in report.source_files["directions"]


def test_source_backed_from_files_rejects_final_physical_file_reuse(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")

    with pytest.raises(ValueError, match="same physical file"):
        RheaDatabase.from_files(release, rdf=release / "raw" / "rhea-directions.tsv")


@pytest.mark.parametrize(
    "role",
    [
        "rdf",
        "directions",
        "relationships",
        "obsolete_reactions",
        "reaction_smiles",
        "sdf",
        "chebi_names",
        "chebi_ph7_3_mapping",
        "xrefs",
        "uniprot_sprot",
        "uniprot_trembl",
    ],
)
def test_source_backed_explicit_roles_reject_missing_and_directory_values(
    tmp_path: Path, role: str
) -> None:
    release = _write_release(tmp_path / f"release-{role}")
    missing = tmp_path / f"missing-{role}"

    with pytest.raises(FileNotFoundError):
        RheaDatabase.from_files(release, **{role: missing})
    with pytest.raises(ValueError, match="not a file"):
        RheaDatabase.from_files(release, **{role: tmp_path})


def test_if_exists_and_failed_build_preserve_destination(
    tmp_path: Path,
) -> None:
    file_names = _write(tmp_path / "names.tsv", "1\t A\n")
    db = RheaDatabase.from_files(chebi_names=file_names)
    path = tmp_path / "rhea.duckdb"
    db.write_duckdb(path)
    original = path.read_bytes()

    with pytest.raises(FileExistsError):
        db.write_duckdb(path)
    assert path.read_bytes() == original

    db.write_duckdb(path, if_exists="replace")
    replacement = path.read_bytes()
    assert replacement

    file_rdf = _write(tmp_path / "rhea.rdf", RDF_FIXTURE)
    file_bad_directions = _write(
        tmp_path / "bad-directions.tsv",
        "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n"
        "20000\t20001\t20002\t20003\n",
    )
    invalid_db = RheaDatabase.from_files(
        rdf=file_rdf,
        directions=file_bad_directions,
    )
    with pytest.raises(ValueError, match="direction semantics disagree"):
        invalid_db.write_duckdb(path, if_exists="replace")
    assert path.read_bytes() == replacement
    assert not list(tmp_path.glob(f".{path.name}.*.tmp*"))


def test_source_hash_is_opt_in(tmp_path: Path) -> None:
    file_names = _write(tmp_path / "names.tsv", "1\t A\n")
    path = tmp_path / "rhea.duckdb"

    RheaDatabase.from_files(chebi_names=file_names).write_duckdb(
        path, include_source_hashes=True
    )

    with duckdb.connect(str(path), read_only=True) as connection:
        row = connection.execute(
            "SELECT sha256 FROM _bioextract.source_file"
        ).fetchone()
    assert row is not None
    digest = row[0]
    assert isinstance(digest, str)
    assert len(digest) == 64


def test_publication_metadata_and_direction_view_follow_shared_contract(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    path = tmp_path / "rhea.duckdb"
    RheaDatabase.from_files(release).write_duckdb(path)

    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        directions = connection.execute(
            """
            SELECT rhea_id, side, directional_role
            FROM reaction_participant_direction
            WHERE participant_id = 'Participant_A'
            ORDER BY rhea_id
            """
        ).fetchall()

    assert metadata["bioextract.resource_name"] == "rhea"
    assert metadata["bioextract.resource_schema_version"] == "rhea-duckdb-v1"
    assert metadata["bioextract.release_version"] == "141"
    assert metadata["bioextract.release_version_source"] == "official_metadata"
    assert "resource_name" not in metadata
    assert directions == [
        (10000, "L", None),
        (10001, "L", "substrate"),
        (10002, "L", "product"),
        (10003, "L", None),
    ]


def test_partial_publication_reports_missing_capabilities(tmp_path: Path) -> None:
    file_names = _write(tmp_path / "names.tsv", "CHEBI:1\t A\n")
    path = tmp_path / "compounds.duckdb"
    RheaDatabase.from_files(chebi_names=file_names).write_duckdb(path)
    database = RheaDatabase.from_duckdb(path)

    with pytest.raises(CapabilityError, match="missing relations"):
        database.select_reactions(["CHEBI:1"], namespace="chebi")
    with pytest.raises(CapabilityError, match="from_duckdb"):
        RheaDatabase.from_files(chebi_names=file_names).select_reactions(
            ["CHEBI:1"], namespace="chebi"
        )
    with pytest.raises(RuntimeError, match="cannot be republished"):
        database.write_duckdb(tmp_path / "copy.duckdb")


def test_from_duckdb_rejects_wrong_identity_and_corrupt_inventory(
    tmp_path: Path,
) -> None:
    wrong = tmp_path / "wrong.duckdb"
    with duckdb.connect(str(wrong)) as connection:
        connection.execute("CREATE TABLE example (value INTEGER)")
    with pytest.raises(ValueError, match="metadata tables"):
        RheaDatabase.from_duckdb(wrong)

    release = _write_release(tmp_path / "release")
    corrupt = tmp_path / "corrupt.duckdb"
    RheaDatabase.from_files(release).write_duckdb(corrupt)
    with duckdb.connect(str(corrupt)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info "
            "SET row_count = row_count + 1 WHERE table_name = 'reaction'"
        )
    with pytest.raises(ValueError, match="row count"):
        RheaDatabase.from_duckdb(corrupt)


def test_metadata_versions_physical_curie_contract_and_native_connection(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    current = tmp_path / "current.duckdb"
    RheaDatabase.from_files(release).write_duckdb(current)

    database = RheaDatabase.from_duckdb(current)
    with database.connect() as connection:
        assert connection.execute(
            "SELECT chebi_id FROM compound WHERE chebi_id IS NOT NULL"
        ).fetchone() == ("CHEBI:1",)
        with pytest.raises(duckdb.Error):
            connection.execute("CREATE TABLE forbidden(value INTEGER)")

    legacy_metadata = tmp_path / "legacy-metadata.duckdb"
    legacy_metadata.write_bytes(current.read_bytes())
    with duckdb.connect(str(legacy_metadata)) as connection:
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
        RheaDatabase.from_duckdb(legacy_metadata)

    unknown = tmp_path / "unknown.duckdb"
    unknown.write_bytes(current.read_bytes())
    with duckdb.connect(str(unknown)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='3' "
            "WHERE key='bioextract.metadata_schema_version'"
        )
    with pytest.raises(ValueError, match="metadata_schema_version"):
        RheaDatabase.from_duckdb(unknown)

    missing_v1_table = tmp_path / "missing-v1-table.duckdb"
    missing_v1_table.write_bytes(current.read_bytes())
    with duckdb.connect(str(missing_v1_table)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="five _bioextract relations"):
        RheaDatabase.from_duckdb(missing_v1_table)

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
            RheaDatabase.from_duckdb(unsupported_profile)

    numeric = tmp_path / "numeric.duckdb"
    numeric.write_bytes(current.read_bytes())
    with duckdb.connect(str(numeric)) as connection:
        connection.execute(
            "ALTER TABLE compound ALTER chebi_id TYPE BIGINT "
            "USING CAST(replace(chebi_id, 'CHEBI:', '') AS BIGINT)"
        )
    with pytest.raises(ValueError, match="VARCHAR CURIE"):
        RheaDatabase.from_duckdb(numeric)
