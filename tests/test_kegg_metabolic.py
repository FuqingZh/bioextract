import io
import tarfile
import zipfile
from pathlib import Path

import duckdb
import pytest

from bioextract._publication import DuckDBWriteResult
from bioextract.kegg import (
    KEGGDatabase,
    KEGGMetabolicCapabilityError,
    KEGGMetabolicNamespace,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def metabolic_files(tmp_path: Path) -> dict[str, Path]:
    compound = _write(
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
    )
    reaction = _write(
        tmp_path / "reaction",
        """ENTRY       R00001                      Reaction
NAME        ATP hydrolysis;
DEFINITION  ATP + H2O to ADP
EQUATION    2 C00001 + n C00002 <=> C00003 + G00001
DBLINKS     Rhea: 12345
RCLASS      RC00001 C00002_C00003
///
""",
    )
    enzyme = _write(
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
    )
    module = _write(
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
    )
    relations = {
        "compound_pubchem": "cpd:C00002\tpubchem:3304\n",
        "compound_reaction": "cpd:C00001\trn:R00001\ncpd:C00002\trn:R00001\n",
        "reaction_enzyme": "rn:R00001\tec:3.6.1.3\n",
        "reaction_ko": "rn:R00001\tko:K00001\n",
        "reaction_module": "rn:R00001\tmd:M00001\n",
        "reaction_pathway": "rn:R00001\tpath:map00010\n",
        "module_pathway": "md:M00001\tpath:map00010\n",
    }
    result = {
        "compound_entries": compound,
        "reaction_entries": reaction,
        "enzyme_entries": enzyme,
        "module_entries": module,
    }
    for name, text in relations.items():
        result[name] = _write(tmp_path / f"{name}.tsv", text)
    return result


def publish(tmp_path: Path) -> tuple[KEGGDatabase, Path]:
    path = tmp_path / "kegg.duckdb"
    source = KEGGDatabase.from_metabolic_files(
        **metabolic_files(tmp_path), release_version="test"
    )
    result = source.write_duckdb(path)
    assert isinstance(result, DuckDBWriteResult)
    assert result.path == path
    return KEGGDatabase.from_duckdb(path), path


def source_from_files(files: dict[str, Path]) -> KEGGDatabase:
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


def test_streaming_parser_and_canonical_relations(tmp_path: Path) -> None:
    db, path = publish(tmp_path)
    with db.connect() as con:
        assert con.execute("SELECT count(*) FROM compound").fetchone() == (2,)
        participants = con.execute(
            "SELECT side, participant_namespace, coefficient_text, coefficient_numeric "
            "FROM reaction_participant ORDER BY side, position"
        ).fetchall()
        assert ("left", "kegg_compound", "n", None) in participants
        assert ("right", "kegg_glycan", "1", 1.0) in participants
        assert con.execute(
            "SELECT external_id FROM compound_cross_reference"
        ).fetchone() == ("CHEBI:15377",)
        assert con.execute(
            "SELECT external_id FROM reaction_cross_reference"
        ).fetchone() == ("RHEA:12345",)
        assert con.execute(
            "SELECT namespace, external_id FROM enzyme_cross_reference"
        ).fetchone() == ("explorenz", "3.6.1.3")
        assert con.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone() == (1,)
    assert path.is_file()


def test_equation_participant_suffix_qualifiers_preserve_coefficients(
    tmp_path: Path,
) -> None:
    reaction = _write(
        tmp_path / "qualified-reaction.keg",
        """ENTRY       R00379                      Reaction
EQUATION    C00039(n+1) + G10477(n) <=> 2 C01330(side 1) + C00001
///
""",
    )
    path = tmp_path / "qualified.duckdb"
    KEGGDatabase.from_metabolic_files(reaction_entries=reaction).write_duckdb(path)

    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute(
            "SELECT participant_id, coefficient_text, coefficient_numeric "
            "FROM reaction_participant ORDER BY side, position"
        ).fetchall() == [
            ("C00039", "1", 1.0),
            ("G10477", "1", 1.0),
            ("C01330", "2", 2.0),
            ("C00001", "1", 1.0),
        ]


def test_read_only_open_selection_and_extractors(tmp_path: Path) -> None:
    db, path = publish(tmp_path)
    selection = db.select_ids(["CHEBI:15377", "CHEBI:999"], namespace="chebi")
    assert selection.extract_matches().select("EntityId").to_series().to_list() == [
        "C00001"
    ]
    assert selection.extract_reactions()["ReactionId"].to_list() == ["R00001"]
    assert selection.extract_participants().height == 3
    assert selection.extract_enzymes()["EcNumber"].to_list() == ["3.6.1.3"]
    assert selection.extract_kos()["KoId"].to_list() == ["K00001"]
    assert selection.extract_modules()["ModuleId"].to_list() == ["M00001"]
    assert selection.extract_pathway_memberships()["PathwayId"].to_list() == [
        "map00010"
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"InputId": "CHEBI:999", "Reason": "not_found"}
    ]
    with db.connect() as con, pytest.raises(duckdb.InvalidInputException):
        con.execute("CREATE TABLE forbidden(i INTEGER)")
    path.write_bytes(path.read_bytes())


def test_grouped_selection_and_exact_module_evaluation(tmp_path: Path) -> None:
    db, _ = publish(tmp_path)
    grouped = db.select_groups(
        {"compound": ["C00002"], "reaction": ["C00001"]},
        namespace="kegg_compound",
    )
    assert set(grouped.extract_reactions()["GroupId"]) == {"compound", "reaction"}
    complete = db.evaluate_modules(["K00001", "K00003"])
    assert complete.to_dicts() == [
        {
            "ModuleId": "M00001",
            "RequiredBlockCount": 2,
            "SatisfiedBlockCount": 2,
            "IsComplete": True,
            "MissingBlockIndexes": [],
        }
    ]
    incomplete = db.evaluate_modules(["K00001"])
    assert incomplete["MissingBlockIndexes"].to_list() == [[2]]


def test_module_double_dash_placeholder_is_optional(tmp_path: Path) -> None:
    module = _write(
        tmp_path / "placeholder-module.keg",
        """ENTRY       M00076                      Pathway module
NAME        Placeholder module
DEFINITION  -- K01136 K01217
REACTION    R00001 C00001 -> C00002

            R00002 C00002 -> C00003
///
""",
    )
    path = tmp_path / "placeholder.duckdb"
    KEGGDatabase.from_metabolic_files(module_entries=module).write_duckdb(path)
    result = KEGGDatabase.from_duckdb(path).evaluate_modules(["K01136", "K01217"])
    assert result.select(
        "RequiredBlockCount", "SatisfiedBlockCount", "IsComplete"
    ).to_dicts() == [
        {
            "RequiredBlockCount": 2,
            "SatisfiedBlockCount": 2,
            "IsComplete": True,
        }
    ]
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute(
            "SELECT reaction_id FROM module_reaction_step ORDER BY position"
        ).fetchall() == [("R00001",), ("R00002",)]


def test_partial_capability_and_publication_validation(tmp_path: Path) -> None:
    files = metabolic_files(tmp_path)
    source = KEGGDatabase.from_metabolic_files(
        reaction_entries=files["reaction_entries"]
    )
    path = tmp_path / "partial.duckdb"
    source.write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    with pytest.raises(KEGGMetabolicCapabilityError, match="namespace 'ko'"):
        db.select_ids(["K00001"], namespace="ko")

    with duckdb.connect(str(path)) as con:
        con.execute(
            "UPDATE _bioextract.metadata SET value='other' "
            "WHERE key='bioextract.resource_name'"
        )
    with pytest.raises(ValueError, match="not a bioextract KEGG"):
        KEGGDatabase.from_duckdb(path)


def test_metadata_inventory_atomic_replace_and_staging_cleanup(
    tmp_path: Path,
) -> None:
    files = metabolic_files(tmp_path)
    source = source_from_files(files)
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


def test_all_metabolic_namespaces_and_curie_cross_references(
    tmp_path: Path,
) -> None:
    db, _ = publish(tmp_path)
    cases: tuple[tuple[KEGGMetabolicNamespace, str], ...] = (
        ("kegg_compound", "C00001"),
        ("chebi", "15377"),
        ("pubchem", "3304"),
        ("kegg_reaction", "R00001"),
        ("rhea", "12345"),
        ("ec", "3.6.1.3"),
        ("ko", "K00001"),
        ("kegg_module", "M00001"),
        ("kegg_pathway", "map00010"),
    )
    for namespace, identifier in cases:
        assert (
            not db.select_ids(
                [identifier],
                namespace=namespace,
            )
            .extract_reactions()
            .is_empty()
        )
    cross_references = db.select_ids(
        ["R00001"], namespace="kegg_reaction"
    ).extract_cross_references()
    assert "RHEA:12345" in cross_references["ExternalId"]
    assert "CHEBI:15377" in cross_references["ExternalId"]


def test_ec_replacement_and_obsolete_reason_precedence(tmp_path: Path) -> None:
    db, _ = publish(tmp_path)
    replacement = db.select_ids(["9.9.9.9"], namespace="ec")
    assert replacement.extract_matches().select("EntityId", "MatchType").to_dicts() == [
        {"EntityId": "3.6.1.3", "MatchType": "replacement"}
    ]
    missing = db.select_ids(["8.8.8.8"], namespace="ec")
    assert missing.extract_unmatched_ids().to_dicts() == [
        {"InputId": "8.8.8.8", "Reason": "not_found"}
    ]


def test_real_enzyme_obsolete_format_and_recursive_replacements(
    tmp_path: Path,
) -> None:
    enzyme = _write(
        tmp_path / "enzyme.keg",
        """ENTRY       EC 1.1.1.1        Obsolete  Enzyme
NAME        Transferred to 1.1.1.2
COMMENT     Transferred entry. Now EC 1.1.1.2
///
ENTRY       EC 1.1.1.2        Obsolete  Enzyme
NAME        Deleted entry
COMMENT     Activity is now covered by EC 1.1.1.3
///
ENTRY       EC 1.1.1.3                  Enzyme
NAME        accepted enzyme;
///
ENTRY       EC 1.1.1.4        Obsolete  Enzyme
NAME        Deleted entry
///
ENTRY       EC 1.1.1.5        Obsolete  Enzyme
NAME        Transferred to 1.1.1.999
///
""",
    )
    path = tmp_path / "enzyme.duckdb"
    KEGGDatabase.from_metabolic_files(enzyme_entries=enzyme).write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)

    assert db.select_ids(["1.1.1.1"], namespace="ec").extract_matches().select(
        "EntityId", "MatchType"
    ).to_dicts() == [{"EntityId": "1.1.1.3", "MatchType": "replacement"}]
    assert db.select_ids(
        ["1.1.1.4"], namespace="ec"
    ).extract_unmatched_ids().to_dicts() == [
        {"InputId": "1.1.1.4", "Reason": "obsolete_excluded"}
    ]
    assert db.select_ids(
        ["1.1.1.5"], namespace="ec"
    ).extract_unmatched_ids().to_dicts() == [
        {"InputId": "1.1.1.5", "Reason": "invalid_canonical_target"}
    ]
    assert db.select_ids(
        ["1.1.1.4"], namespace="ec", include_obsolete=True
    ).extract_matches().select("EntityId", "MatchType").to_dicts() == [
        {"EntityId": "1.1.1.4", "MatchType": "exact"}
    ]
    with db.connect() as connection:
        assert connection.execute(
            "SELECT ec_number, status FROM enzyme ORDER BY ec_number"
        ).fetchall() == [
            ("1.1.1.1", "transferred"),
            ("1.1.1.2", "transferred"),
            ("1.1.1.3", "active"),
            ("1.1.1.4", "deleted"),
            ("1.1.1.5", "transferred"),
        ]
        assert connection.execute(
            "SELECT ec_number, replacement_ec_number "
            "FROM enzyme_replacement ORDER BY ec_number"
        ).fetchall() == [
            ("1.1.1.1", "1.1.1.2"),
            ("1.1.1.2", "1.1.1.3"),
        ]
        assert connection.execute(
            "SELECT issue_code, referenced_identifier FROM _bioextract.validation_issue"
        ).fetchall() == [("foreign_key_violation", "1.1.1.999")]


def test_release_discovery_uses_exact_layout(tmp_path: Path) -> None:
    raw = tmp_path / "2026-07" / "raw"
    files = metabolic_files(tmp_path)
    for family in ("compound", "reaction", "enzyme", "module"):
        entries = raw / family / "entries"
        entries.mkdir(parents=True)
        (raw / family / "list.tsv").write_text(
            f"{family[:1].upper()}00001\tentry\n", encoding="utf-8"
        )
        (entries / "000001.keg").write_text(
            files[f"{family}_entries"].read_text(encoding="utf-8"),
            encoding="utf-8",
        )
    for relation in (
        "compound_pubchem",
        "compound_reaction",
        "reaction_enzyme",
        "reaction_ko",
        "reaction_module",
        "reaction_pathway",
        "module_pathway",
    ):
        (raw / f"{relation}.tsv").write_text(
            files[relation].read_text(encoding="utf-8"), encoding="utf-8"
        )
    (raw / "reaction" / "ignore.tsv").write_text("noise", encoding="utf-8")
    db = KEGGDatabase.from_metabolic_release(raw.parent)
    snapshot = db.snapshot.metabolic
    assert snapshot is not None
    assert len(snapshot.sources["reaction_entries"]) == 1
    assert snapshot.sources["reaction_list"] == (raw / "reaction" / "list.tsv",)


def test_release_archive_accepts_an_extra_top_level_directory(tmp_path: Path) -> None:
    files = metabolic_files(tmp_path)
    raw = tmp_path / "archive-root" / "wrapper" / "2026-07" / "raw"
    list_ids = {
        "compound": ("C00001", "C00002"),
        "reaction": ("R00001",),
        "enzyme": ("3.6.1.3", "9.9.9.9"),
        "module": ("M00001",),
    }
    for family, identifiers in list_ids.items():
        entries = raw / family / "entries"
        entries.mkdir(parents=True)
        (entries / f"{family}.keg").write_bytes(files[f"{family}_entries"].read_bytes())
        (raw / family / "list.tsv").write_text(
            "".join(f"{identifier}\tentry\n" for identifier in identifiers),
            encoding="utf-8",
        )
    for relation in (
        "compound_pubchem",
        "compound_reaction",
        "reaction_enzyme",
        "reaction_ko",
        "reaction_module",
        "reaction_pathway",
        "module_pathway",
    ):
        (raw / f"{relation}.tsv").write_bytes(files[relation].read_bytes())

    archive = tmp_path / "release.zip"
    with zipfile.ZipFile(archive, "w") as output:
        for source in sorted((tmp_path / "archive-root").rglob("*")):
            if source.is_file():
                output.write(source, source.relative_to(tmp_path / "archive-root"))
    path = tmp_path / "release.duckdb"
    KEGGDatabase.from_metabolic_release(archive).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute("SELECT count(*) FROM compound").fetchone() == (2,)


def test_existing_kegg_modes_reject_metabolic_only_operations(tmp_path: Path) -> None:
    file_brite = _write(tmp_path / "brite.json", '{"name":"x","children":[]}')
    db = KEGGDatabase.from_brite_json(file_brite)
    with pytest.raises(KEGGMetabolicCapabilityError):
        db.connect()


def test_capability_metadata_matches_actual_inventory(tmp_path: Path) -> None:
    _, path = publish(tmp_path)
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


def test_metadata_v3_requires_validation_issue_table(tmp_path: Path) -> None:
    _, path = publish(tmp_path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="validation_issue"):
        KEGGDatabase.from_duckdb(path)


def test_relation_only_inputs_preserve_rows_and_enable_namespace(
    tmp_path: Path,
) -> None:
    relation = _write(tmp_path / "reaction_enzyme.tsv", "rn:R00001\tec:3.6.1.3\n")
    path = tmp_path / "relation-only.duckdb"
    KEGGDatabase.from_metabolic_files(reaction_enzyme=relation).write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    selection = db.select_ids(["3.6.1.3"], namespace="ec")
    assert selection.extract_matches().to_dicts() == [
        {
            "InputId": "3.6.1.3",
            "InputNamespace": "ec",
            "EntityType": "reaction",
            "EntityId": "R00001",
            "MatchType": "exact",
        }
    ]


def test_reaction_only_compound_selection_uses_participants(
    tmp_path: Path,
) -> None:
    reaction = _write(
        tmp_path / "reaction.keg",
        "ENTRY       R00001                      Reaction\n"
        "EQUATION    C00001 <=> C00002\n///\n",
    )
    path = tmp_path / "reaction-only.duckdb"
    KEGGDatabase.from_metabolic_files(reaction_entries=reaction).write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    assert db.select_ids(["C00001"], namespace="kegg_compound").extract_reactions()[
        "ReactionId"
    ].to_list() == ["R00001"]


def test_cross_reference_namespace_must_be_present_in_rows(tmp_path: Path) -> None:
    compound = _write(
        tmp_path / "compound.keg",
        "ENTRY       C00001                      Compound\n"
        "DBLINKS     PubChem: 123\n///\n",
    )
    path = tmp_path / "compound.duckdb"
    KEGGDatabase.from_metabolic_files(compound_entries=compound).write_duckdb(path)
    db = KEGGDatabase.from_duckdb(path)
    with pytest.raises(KEGGMetabolicCapabilityError, match="namespace 'chebi'"):
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


def test_module_references_are_recursive_and_cycles_fail(tmp_path: Path) -> None:
    modules = _write(
        tmp_path / "modules.keg",
        "ENTRY       M00001                      Pathway module\n"
        "DEFINITION  K00001 M00002\n///\n"
        "ENTRY       M00002                      Pathway module\n"
        "DEFINITION  K00002\n///\n",
    )
    path = tmp_path / "modules.duckdb"
    KEGGDatabase.from_metabolic_files(module_entries=modules).write_duckdb(path)
    result = KEGGDatabase.from_duckdb(path).evaluate_modules(["K00001", "K00002"])
    assert result.filter(result["ModuleId"] == "M00001")["IsComplete"].item()

    cyclic = _write(
        tmp_path / "cyclic.keg",
        "ENTRY       M00001                      Pathway module\n"
        "DEFINITION  M00002\n///\n"
        "ENTRY       M00002                      Pathway module\n"
        "DEFINITION  M00001\n///\n",
    )
    cycle_path = tmp_path / "cycle.duckdb"
    KEGGDatabase.from_metabolic_files(module_entries=cyclic).write_duckdb(cycle_path)
    with pytest.raises(ValueError, match="Cyclic KEGG module"):
        KEGGDatabase.from_duckdb(cycle_path).evaluate_modules([])
