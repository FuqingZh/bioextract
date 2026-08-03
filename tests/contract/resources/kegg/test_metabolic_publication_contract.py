import io
import tarfile
import zipfile
from pathlib import Path

import duckdb
import pytest

from bioextract.errors import CapabilityError
from bioextract.kegg import KEGGDatabase


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _metabolic_files(tmp_path: Path) -> dict[str, Path]:
    result = {
        "compound_entries": _write(
            tmp_path / "compound",
            """ENTRY       C00001                      Compound
NAME        Water;
FORMULA     H2O
EXACT_MASS  18.0106
MOL_WEIGHT  18.015
DBLINKS     ChEBI: 15377
///
ENTRY       C00002                      Compound
NAME        ATP;
DBLINKS     PubChem: 3304
///
""",
        ),
        "reaction_entries": _write(
            tmp_path / "reaction",
            """ENTRY       R00001                      Reaction
NAME        ATP hydrolysis;
DEFINITION  ATP + H2O to ADP
EQUATION    2 C00001 + n C00002 <=> C00003 + G00001
DBLINKS     Rhea: 12345
RCLASS      RC00001 C00002_C00003
///
""",
        ),
        "enzyme_entries": _write(
            tmp_path / "enzyme",
            """ENTRY       EC 3.6.1.3                  Enzyme
NAME        adenosinetriphosphatase;
CLASS       Hydrolases
ORTHOLOGY   K00001 ATPase
DBLINKS     ExplorEnz: 3.6.1.3
///
ENTRY       EC 9.9.9.9                  Enzyme
DEFINITION  Transferred entry
TRANSFER    3.6.1.3
///
""",
        ),
        "module_entries": _write(
            tmp_path / "module",
            """ENTRY       M00001                      Pathway module
NAME        Test module
CLASS       Metabolism
DEFINITION  K00001 (K00002,K00003) -K00004
REACTION    R00001 C00002 -> C00003
COMPOUND    C00002 ATP
DIAGRAM     M00001
///
""",
        ),
    }
    relations = {
        "compound_pubchem": "cpd:C00002\tpubchem:3304\n",
        "compound_reaction": "cpd:C00001\trn:R00001\ncpd:C00002\trn:R00001\n",
        "reaction_enzyme": "rn:R00001\tec:3.6.1.3\n",
        "reaction_ko": "rn:R00001\tko:K00001\n",
        "reaction_module": "rn:R00001\tmd:M00001\n",
        "reaction_pathway": "rn:R00001\tpath:map00010\n",
        "module_pathway": "md:M00001\tpath:map00010\n",
    }
    result.update(
        {
            name: _write(tmp_path / f"{name}.tsv", text)
            for name, text in relations.items()
        }
    )
    return result


def _source_from_files(files: dict[str, Path]) -> KEGGDatabase:
    return KEGGDatabase.from_metabolic_files(
        compound_entries=files["compound_entries"],
        reaction_entries=files["reaction_entries"],
        enzyme_entries=files["enzyme_entries"],
        module_entries=files["module_entries"],
        compound_pubchem=files["compound_pubchem"],
        compound_reaction=files["compound_reaction"],
        reaction_enzyme=files["reaction_enzyme"],
        reaction_ko=files["reaction_ko"],
        reaction_module=files["reaction_module"],
        reaction_pathway=files["reaction_pathway"],
        module_pathway=files["module_pathway"],
    )


def _publish(tmp_path: Path) -> Path:
    path = tmp_path / "kegg.duckdb"
    KEGGDatabase.from_metabolic_files(
        **_metabolic_files(tmp_path),
        release_version="test",
    ).write_duckdb(path)
    return path


def test_partial_capability_and_publication_validation(tmp_path: Path) -> None:
    files = _metabolic_files(tmp_path)
    source = KEGGDatabase.from_metabolic_files(
        reaction_entries=files["reaction_entries"]
    )
    path = tmp_path / "partial.duckdb"
    source.write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    with pytest.raises(CapabilityError, match="namespace 'ko'"):
        db.select_ids(["K00001"], namespace="ko")

    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='other' "
            "WHERE key='bioextract.resource_name'"
        )
    with pytest.raises(ValueError, match="not a bioextract KEGG"):
        KEGGDatabase.from_duckdb(path)


def test_metadata_inventory_atomic_replace_and_staging_cleanup(
    tmp_path: Path,
) -> None:
    files = _metabolic_files(tmp_path)
    source = _source_from_files(files)
    path = tmp_path / "kegg.duckdb"
    source.write_duckdb(path)
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        source.write_duckdb(path)
    assert path.read_bytes() == original
    assert not list(tmp_path.glob(f".{path.name}.*.duckdb"))

    with duckdb.connect(str(path), read_only=True) as connection:
        metadata_tables = {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema='_bioextract'"
            ).fetchall()
        }
        assert metadata_tables == {
            "metadata",
            "source_file",
            "table_info",
            "column_mapping",
            "validation_issue",
        }
        for table, count in connection.execute(
            "SELECT table_name, row_count FROM _bioextract.table_info"
        ).fetchall():
            assert connection.execute(f'SELECT count(*) FROM "{table}"').fetchone() == (
                count,
            )


def test_existing_kegg_modes_reject_metabolic_only_operations(tmp_path: Path) -> None:
    file_brite = _write(tmp_path / "brite.json", '{"name":"x","children":[]}')
    db = KEGGDatabase.from_brite_json(file_brite)
    with pytest.raises(CapabilityError):
        db.connect()


def test_capability_metadata_matches_actual_inventory(tmp_path: Path) -> None:
    path = _publish(tmp_path)
    with duckdb.connect(str(path)) as connection:
        capabilities = connection.execute(
            "SELECT value FROM _bioextract.metadata WHERE key='bioextract.capabilities'"
        ).fetchone()
        assert capabilities is not None
        assert "enzyme_entries" in capabilities[0].split(",")
        connection.execute(
            "UPDATE _bioextract.metadata SET value=value || ',unknown_role' "
            "WHERE key='bioextract.capabilities'"
        )
    with pytest.raises(ValueError, match="unknown capabilities"):
        KEGGDatabase.from_duckdb(path)


def test_metadata_v1_requires_validation_issue_table(tmp_path: Path) -> None:
    path = _publish(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="validation_issue"):
        KEGGDatabase.from_duckdb(path)


@pytest.mark.parametrize("profile", ["", "unknown-profile"])
def test_metadata_v1_requires_supported_source_schema_profile(
    tmp_path: Path, profile: str
) -> None:
    path = _publish(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value=? "
            "WHERE key='bioextract.source_schema_profile'",
            [profile],
        )
    with pytest.raises(ValueError, match="source schema profile"):
        KEGGDatabase.from_duckdb(path)


def test_cross_reference_namespace_must_be_present_in_rows(tmp_path: Path) -> None:
    compound = _write(
        tmp_path / "compound.keg",
        "ENTRY       C00001                      Compound\n"
        "DBLINKS     PubChem: 123\n///\n",
    )
    path = tmp_path / "compound.duckdb"
    KEGGDatabase.from_metabolic_files(compound_entries=compound).write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    with pytest.raises(CapabilityError, match="namespace 'chebi'"):
        db.select_ids(["CHEBI:15377"], namespace="chebi")


def test_entity_presence_controls_orphan_validation(tmp_path: Path) -> None:
    compound = _write(
        tmp_path / "compound.keg",
        "ENTRY       C00001                      Compound\n///\n",
    )
    relation = _write(
        tmp_path / "compound_reaction.tsv",
        "cpd:C00001\trn:R00001\ncpd:C99999\trn:R99999\n",
    )
    path = tmp_path / "validated.duckdb"
    KEGGDatabase.from_metabolic_files(
        compound_entries=compound, compound_reaction=relation
    ).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute(
            "SELECT compound_id, reaction_id FROM compound_reaction"
        ).fetchall() == [("C00001", "R00001")]
        assert connection.execute(
            "SELECT issue_code FROM _bioextract.validation_issue"
        ).fetchall() == [("foreign_key_violation",)]


def test_missing_compound_participant_is_skipped_only_with_compound_entities(
    tmp_path: Path,
) -> None:
    reaction = _write(
        tmp_path / "reaction.keg",
        "ENTRY       R00001                      Reaction\n"
        "EQUATION    C00001 <=> C99999\n///\n",
    )
    compound = _write(
        tmp_path / "compound.keg",
        "ENTRY       C00001                      Compound\n///\n",
    )
    partial = tmp_path / "partial.duckdb"
    KEGGDatabase.from_metabolic_files(reaction_entries=reaction).write_duckdb(partial)
    with KEGGDatabase.from_duckdb(partial).connect() as connection:
        assert connection.execute(
            "SELECT participant_id FROM reaction_participant ORDER BY position"
        ).fetchall() == [("C00001",), ("C99999",)]

    validated = tmp_path / "with-compounds.duckdb"
    KEGGDatabase.from_metabolic_files(
        compound_entries=compound, reaction_entries=reaction
    ).write_duckdb(validated)
    with KEGGDatabase.from_duckdb(validated).connect() as connection:
        assert connection.execute(
            "SELECT participant_id FROM reaction_participant"
        ).fetchall() == [("C00001",)]
        assert connection.execute(
            "SELECT issue_code, identifier_value, referenced_identifier "
            "FROM _bioextract.validation_issue"
        ).fetchall() == [("foreign_key_violation", "R00001", "C99999")]


@pytest.mark.parametrize("kind", ["zip", "tar"])
def test_entry_archives_are_streamed_safely(tmp_path: Path, kind: str) -> None:
    content = b"ENTRY       C00001                      Compound\n///\n"
    archive = tmp_path / f"compound.{kind}"
    if kind == "zip":
        with zipfile.ZipFile(archive, "w") as output:
            output.writestr("wrapper/entries/compound.keg", content)
    else:
        with tarfile.open(archive, "w") as output:
            info = tarfile.TarInfo("wrapper/entries/compound.keg")
            info.size = len(content)
            output.addfile(info, io.BytesIO(content))
    path = tmp_path / f"{kind}.duckdb"
    KEGGDatabase.from_metabolic_files(compound_entries=archive).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute("SELECT compound_id FROM compound").fetchall() == [
            ("C00001",)
        ]
        source = connection.execute(
            "SELECT display_path, media_type FROM _bioextract.source_file"
        ).fetchone()
        assert source is not None
        assert source[0] == str(archive)
        assert source[1] in {"application/zip", "application/x-tar"}


@pytest.mark.parametrize("kind", ["zip", "tar"])
def test_entry_archive_rejects_path_traversal(tmp_path: Path, kind: str) -> None:
    content = b"ENTRY       C00001                      Compound\n///\n"
    archive = tmp_path / f"unsafe.{kind}"
    if kind == "zip":
        with zipfile.ZipFile(archive, "w") as output:
            output.writestr("../escape.keg", content)
    else:
        with tarfile.open(archive, "w") as output:
            info = tarfile.TarInfo("../escape.keg")
            info.size = len(content)
            output.addfile(info, io.BytesIO(content))
    with pytest.raises(ValueError, match="Unsafe path"):
        KEGGDatabase.from_metabolic_files(compound_entries=archive).write_duckdb(
            tmp_path / "unsafe.duckdb"
        )
    assert not (tmp_path / "escape.keg").exists()
