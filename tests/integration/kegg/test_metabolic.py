import shutil
import tarfile
import zipfile
from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract._publication import DuckDBWriteResult
from bioextract.kegg import KEGGDatabase
from bioextract.kegg.metabolic.core import (
    KEGGMetabolicNamespace,
    KEGGMetabolicSelection,
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


def metabolic_source(
    files: dict[str, Path],
    *,
    source: Path | None = None,
    release_version: str | None = None,
) -> KEGGDatabase:
    return KEGGDatabase.from_metabolic_files(
        source,
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
        release_version=release_version,
    )


def publish(tmp_path: Path) -> tuple[KEGGDatabase, Path]:
    path = tmp_path / "kegg.duckdb"
    source = metabolic_source(
        metabolic_files(tmp_path),
        release_version="test",
    )
    result = source.write_duckdb(path)
    assert isinstance(result, DuckDBWriteResult)
    assert result.path == path
    return KEGGDatabase.from_duckdb(path), path


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


def test_read_only_open_selection_and_extractors(tmp_path: Path) -> None:
    db, path = publish(tmp_path)
    selection = db.select_ids(["CHEBI:15377", "CHEBI:999"], namespace="chebi")
    assert selection.extract_matches().select("entity_id").to_series().to_list() == [
        "C00001"
    ]
    assert selection.extract_reactions()["reaction_id"].to_list() == ["R00001"]
    assert selection.extract_participants().height == 3
    assert selection.extract_enzymes()["ec_number"].to_list() == ["3.6.1.3"]
    assert selection.extract_kos()["ko_id"].to_list() == ["K00001"]
    assert selection.extract_modules()["module_id"].to_list() == ["M00001"]
    assert selection.extract_pathway_memberships()["pathway_id"].to_list() == [
        "map00010"
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"input_id": "CHEBI:999", "reason": "not_found"}
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
    assert set(grouped.extract_reactions()["group_id"]) == {"compound", "reaction"}
    complete = db.evaluate_modules(["K00001", "K00003"])
    assert complete.to_dicts() == [
        {
            "module_id": "M00001",
            "required_block_count": 2,
            "satisfied_block_count": 2,
            "is_complete": True,
            "missing_block_indexes": [],
        }
    ]
    incomplete = db.evaluate_modules(["K00001"])
    assert incomplete["missing_block_indexes"].to_list() == [[2]]


def test_grouped_selection_resolves_unique_ids_once_then_expands(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db, _ = publish(tmp_path)
    calls: list[tuple[str, ...]] = []
    original = KEGGMetabolicSelection._query_unique_matches  # pyright: ignore[reportPrivateUsage]

    def tracked_query(
        selection: KEGGMetabolicSelection,
        connection: duckdb.DuckDBPyConnection,
    ) -> pl.DataFrame:
        calls.append(selection.input_ids)
        return original(selection, connection)

    monkeypatch.setattr(
        KEGGMetabolicSelection,
        "_query_unique_matches",
        tracked_query,
    )
    selection = db.select_groups(
        {
            " first ": ["cpd:C00001", "C00001", ""],
            "second": ["C00001", "C99999", "C99999"],
            "empty": [" "],
        },
        namespace="kegg_compound",
    )

    assert selection.input_ids == ("C00001", "C99999")
    assert selection.group_ids == ("empty", "first", "second")
    assert selection.group_membership == (
        ("first", "C00001"),
        ("second", "C00001"),
        ("second", "C99999"),
    )
    assert selection.extract_matches().select("group_id", "input_id").to_dicts() == [
        {"group_id": "first", "input_id": "C00001"},
        {"group_id": "second", "input_id": "C00001"},
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [
        {"group_id": "second", "input_id": "C99999", "reason": "not_found"}
    ]
    assert set(selection.extract_reactions()["group_id"]) == {"first", "second"}

    empty_selection = db.select_groups(
        {"empty": [" "]},
        namespace="kegg_compound",
    )
    assert empty_selection.is_grouped
    assert empty_selection.group_ids == ("empty",)
    assert empty_selection.extract_matches().columns[0] == "group_id"
    assert empty_selection.extract_unmatched_ids().columns[0] == "group_id"
    assert calls == [("C00001", "C99999")]


def test_grouped_selection_rejects_colliding_normalized_group_ids(
    tmp_path: Path,
) -> None:
    db, _ = publish(tmp_path)
    with pytest.raises(ValueError, match="unique after normalization"):
        db.select_groups(
            {"group": ["C00001"], " group ": ["C00002"]},
            namespace="kegg_compound",
        )


def test_module_placeholder_round_trips_through_evaluation_and_reaction_steps(
    tmp_path: Path,
) -> None:
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
        "required_block_count", "satisfied_block_count", "is_complete"
    ).to_dicts() == [
        {
            "required_block_count": 2,
            "satisfied_block_count": 2,
            "is_complete": True,
        }
    ]
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute(
            "SELECT reaction_id FROM module_reaction_step ORDER BY position"
        ).fetchall() == [("R00001",), ("R00002",)]


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
    assert "RHEA:12345" in cross_references["external_id"]
    assert "CHEBI:15377" in cross_references["external_id"]


def test_ec_replacement_and_obsolete_reason_precedence(tmp_path: Path) -> None:
    db, _ = publish(tmp_path)
    replacement = db.select_ids(["9.9.9.9"], namespace="ec")
    assert replacement.extract_matches().select(
        "entity_id", "match_type"
    ).to_dicts() == [{"entity_id": "3.6.1.3", "match_type": "replacement"}]
    missing = db.select_ids(["8.8.8.8"], namespace="ec")
    assert missing.extract_unmatched_ids().to_dicts() == [
        {"input_id": "8.8.8.8", "reason": "not_found"}
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
        "entity_id", "match_type"
    ).to_dicts() == [{"entity_id": "1.1.1.3", "match_type": "replacement"}]
    assert db.select_ids(
        ["1.1.1.4"], namespace="ec"
    ).extract_unmatched_ids().to_dicts() == [
        {"input_id": "1.1.1.4", "reason": "obsolete_excluded"}
    ]
    assert db.select_ids(
        ["1.1.1.5"], namespace="ec"
    ).extract_unmatched_ids().to_dicts() == [
        {"input_id": "1.1.1.5", "reason": "invalid_canonical_target"}
    ]
    assert db.select_ids(
        ["1.1.1.4"], namespace="ec", include_obsolete=True
    ).extract_matches().select("entity_id", "match_type").to_dicts() == [
        {"entity_id": "1.1.1.4", "match_type": "exact"}
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


def metabolic_release(tmp_path: Path, *, directory_name: str = "2026-07") -> Path:
    release = tmp_path / directory_name
    raw = release / "raw"
    files = metabolic_files(tmp_path)
    list_ids = {
        "compound": ("C00001", "C00002"),
        "reaction": ("R00001",),
        "enzyme": ("3.6.1.3", "9.9.9.9"),
        "module": ("M00001",),
    }
    for family in ("compound", "reaction", "enzyme", "module"):
        entries = raw / family / "entries"
        entries.mkdir(parents=True)
        (raw / family / "list.tsv").write_text(
            "".join(f"{identifier}\tentry\n" for identifier in list_ids[family]),
            encoding="utf-8",
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
    return release


def _archive_release(
    release: Path,
    archive: Path,
    *,
    kind: str,
    prefix: Path = Path("wrapper") / "2026-07",
) -> Path:
    if kind == "zip":
        with zipfile.ZipFile(archive, "w") as output:
            for source in sorted(release.rglob("*")):
                if source.is_file():
                    output.write(source, prefix / source.relative_to(release))
    else:
        with tarfile.open(archive, "w") as output:
            for source in sorted(release.rglob("*")):
                if source.is_file():
                    output.add(
                        source,
                        arcname=str(prefix / source.relative_to(release)),
                    )
    return archive


@pytest.mark.parametrize("layout", ["release", "raw", "nested"])
def test_release_discovery_uses_one_exact_layout(
    tmp_path: Path,
    layout: str,
) -> None:
    release = metabolic_release(tmp_path)
    raw = release / "raw"
    source = release
    if layout == "raw":
        source = raw
    elif layout == "nested":
        source = tmp_path / "nested-source"
        raw = source / "wrapper" / "raw"
        shutil.copytree(release / "raw", raw)

    db = KEGGDatabase.from_metabolic_files(source)
    snapshot = db.snapshot.metabolic
    assert snapshot is not None
    assert snapshot.complete_release is True
    assert len(snapshot.sources["reaction_entries"]) == 1
    assert snapshot.sources["reaction_list"] == (raw / "reaction" / "list.tsv",)


def test_release_directory_name_does_not_create_release_metadata(
    tmp_path: Path,
) -> None:
    release = metabolic_release(tmp_path, directory_name="2099-12")
    path = tmp_path / "unknown-release.duckdb"
    KEGGDatabase.from_metabolic_files(release).write_duckdb(path)

    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    assert "bioextract.release_version" not in metadata
    assert "bioextract.release_version_source" not in metadata


def test_caller_release_version_is_recorded_with_caller_source(tmp_path: Path) -> None:
    release = metabolic_release(tmp_path, directory_name="arbitrary-layout")
    path = tmp_path / "caller-release.duckdb"
    KEGGDatabase.from_metabolic_files(release, release_version="2026-07").write_duckdb(
        path
    )

    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    assert metadata["bioextract.release_version"] == "2026-07"
    assert metadata["bioextract.release_version_source"] == "caller"


@pytest.mark.parametrize("kind", ["zip", "tar"])
def test_release_archive_accepts_nested_layout(
    tmp_path: Path,
    kind: str,
) -> None:
    release = metabolic_release(tmp_path)
    archive = _archive_release(release, tmp_path / f"release.{kind}", kind=kind)
    path = tmp_path / "release.duckdb"
    KEGGDatabase.from_metabolic_files(archive).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute("SELECT count(*) FROM compound").fetchone() == (2,)


def test_source_overlays_replace_whole_roles_and_final_provenance(
    tmp_path: Path,
) -> None:
    release = metabolic_release(tmp_path)
    raw = release / "raw"
    original_entries = raw / "reaction" / "entries" / "000001.keg"
    original_relation = raw / "reaction_ko.tsv"
    overlay_entries = tmp_path / "overlay-reaction"
    overlay_entries.mkdir()
    overlay_entry = overlay_entries / "replacement.keg"
    overlay_entry.write_bytes(original_entries.read_bytes())
    overlay_relation = tmp_path / "overlay-reaction-ko.tsv"
    overlay_relation.write_bytes(original_relation.read_bytes())
    original_relation.unlink()

    database = KEGGDatabase.from_metabolic_files(
        release,
        reaction_entries=overlay_entries,
        reaction_ko=overlay_relation,
    )
    snapshot = database.snapshot.metabolic
    assert snapshot is not None
    assert snapshot.complete_release is True
    assert snapshot.sources["reaction_entries"] == (overlay_entry,)
    assert snapshot.sources["reaction_ko"] == (overlay_relation,)

    path = tmp_path / "overlay.duckdb"
    database.write_duckdb(path)
    with duckdb.connect(str(path), read_only=True) as connection:
        display_paths = {
            row[0]
            for row in connection.execute(
                "SELECT display_path FROM _bioextract.source_file"
            ).fetchall()
        }
    assert str(overlay_entry) in display_paths
    assert str(overlay_relation) in display_paths
    assert str(original_entries) not in display_paths
    assert str(original_relation) not in display_paths


def test_source_rejects_a_missing_unoverridden_required_role(tmp_path: Path) -> None:
    release = metabolic_release(tmp_path)
    (release / "raw" / "reaction_ko.tsv").unlink()

    with pytest.raises(ValueError, match="missing roles.*reaction_ko"):
        KEGGDatabase.from_metabolic_files(release)


@pytest.mark.parametrize(
    ("role", "value_kind"),
    [
        ("compound_list", "directory"),
        ("reaction_entries", "directory"),
        ("reaction_entries", "sequence"),
        ("reaction_ko", "directory"),
    ],
)
def test_explicit_non_none_role_rejects_an_empty_fileset(
    tmp_path: Path,
    role: str,
    value_kind: str,
) -> None:
    release = metabolic_release(tmp_path)
    empty_directory = tmp_path / f"empty-{role}"
    if value_kind == "directory":
        empty_directory.mkdir()

    with pytest.raises(ValueError, match=rf"role '{role}'.*at least one file"):
        if role == "compound_list":
            KEGGDatabase.from_metabolic_files(
                release,
                compound_list=empty_directory,
            )
        elif role == "reaction_entries" and value_kind == "sequence":
            KEGGDatabase.from_metabolic_files(release, reaction_entries=[])
        elif role == "reaction_entries":
            KEGGDatabase.from_metabolic_files(
                release,
                reaction_entries=empty_directory,
            )
        else:
            KEGGDatabase.from_metabolic_files(
                release,
                reaction_ko=empty_directory,
            )


@pytest.mark.parametrize("role", ["compound_list", "reaction_ko"])
def test_scalar_role_accepts_a_nonempty_directory(
    tmp_path: Path,
    role: str,
) -> None:
    directory = tmp_path / role
    directory.mkdir()
    source = _write(directory / "source.tsv", "value\n")

    database = (
        KEGGDatabase.from_metabolic_files(compound_list=directory)
        if role == "compound_list"
        else KEGGDatabase.from_metabolic_files(reaction_ko=directory)
    )
    snapshot = database.snapshot.metabolic
    assert snapshot is not None
    assert snapshot.sources[role] == (source,)


def test_entry_collection_order_is_deterministic(tmp_path: Path) -> None:
    directory = tmp_path / "compound-entries"
    directory.mkdir()
    second = _write(directory / "b.keg", "ENTRY       C00002\n///\n")
    first = _write(directory / "a.keg", "ENTRY       C00001\n///\n")

    database = KEGGDatabase.from_metabolic_files(compound_entries=directory)
    snapshot = database.snapshot.metabolic
    assert snapshot is not None
    assert snapshot.sources["compound_entries"] == (first, second)

    sequence = KEGGDatabase.from_metabolic_files(
        compound_entries=[second, first]
    ).snapshot.metabolic
    assert sequence is not None
    assert sequence.sources["compound_entries"] == (first, second)


def test_complete_release_allows_absent_optional_lists(tmp_path: Path) -> None:
    release = metabolic_release(tmp_path)
    for file_list in release.glob("raw/*/list.tsv"):
        file_list.unlink()

    path = tmp_path / "without-lists.duckdb"
    KEGGDatabase.from_metabolic_files(release).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute("SELECT count(*) FROM compound").fetchone() == (2,)


def test_present_optional_list_still_requires_entry_parity(tmp_path: Path) -> None:
    release = metabolic_release(tmp_path)
    (release / "raw" / "compound" / "list.tsv").write_text(
        "C99999\tentry\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="compound list/entry mismatch"):
        KEGGDatabase.from_metabolic_files(release).write_duckdb(
            tmp_path / "mismatch.duckdb"
        )


def test_source_layout_is_required_even_when_all_required_roles_are_overlaid(
    tmp_path: Path,
) -> None:
    files = metabolic_files(tmp_path)
    empty_source = tmp_path / "empty-source"
    empty_source.mkdir()

    with pytest.raises(ValueError, match="contains no release layout"):
        metabolic_source(files, source=empty_source)


def test_source_rejects_multiple_plausible_layouts(tmp_path: Path) -> None:
    files = metabolic_files(tmp_path)
    source = tmp_path / "ambiguous-source"
    for prefix in ("first", "second"):
        directory = source / prefix / "raw" / "compound"
        directory.mkdir(parents=True)
        (directory / "list.tsv").write_text("C00001\tentry\n", encoding="utf-8")

    with pytest.raises(ValueError, match="multiple release layouts"):
        metabolic_source(files, source=source)


@pytest.mark.parametrize("kind", ["zip", "tar"])
def test_archive_rejects_zero_or_multiple_plausible_layouts(
    tmp_path: Path,
    kind: str,
) -> None:
    files = metabolic_files(tmp_path)
    empty_source = tmp_path / "empty-archive-source"
    empty_source.mkdir()
    _write(empty_source / "README.txt", "not a KEGG release\n")
    empty_archive = _archive_release(
        empty_source,
        tmp_path / f"empty.{kind}",
        kind=kind,
        prefix=Path("."),
    )
    with pytest.raises(ValueError, match="contains no release layout"):
        metabolic_source(files, source=empty_archive)

    ambiguous_source = tmp_path / "ambiguous-archive-source"
    for prefix in ("first", "second"):
        directory = ambiguous_source / prefix / "raw" / "compound"
        directory.mkdir(parents=True)
        (directory / "list.tsv").write_text("C00001\tentry\n", encoding="utf-8")
    ambiguous_archive = _archive_release(
        ambiguous_source,
        tmp_path / f"ambiguous.{kind}",
        kind=kind,
        prefix=Path("."),
    )
    with pytest.raises(ValueError, match="multiple release layouts"):
        metabolic_source(files, source=ambiguous_archive)


def test_archive_rejects_duplicate_normalized_members(tmp_path: Path) -> None:
    files = metabolic_files(tmp_path)
    archive = tmp_path / "duplicate-members.zip"
    with zipfile.ZipFile(archive, "w") as output:
        output.writestr("raw/compound/list.tsv", "C00001\tentry\n")
        output.writestr("./raw/compound/list.tsv", "C00001\tentry\n")

    with pytest.raises(ValueError, match="duplicate members"):
        metabolic_source(files, source=archive)


def test_final_inventory_rejects_physical_file_reuse(tmp_path: Path) -> None:
    source = _write(tmp_path / "relation.tsv", "rn:R00001\tko:K00001\n")

    with pytest.raises(ValueError, match="same physical file"):
        KEGGDatabase.from_metabolic_files(
            reaction_ko=source,
            reaction_module=source,
        )


def test_release_archive_cannot_also_be_an_explicit_role(tmp_path: Path) -> None:
    release = metabolic_release(tmp_path)
    archive = _archive_release(release, tmp_path / "release.zip", kind="zip")

    with pytest.raises(ValueError, match="same physical file"):
        KEGGDatabase.from_metabolic_files(archive, reaction_entries=archive)


@pytest.mark.parametrize("kind", ["zip", "tar"])
def test_release_archive_overlay_satisfies_missing_role_and_records_provenance(
    tmp_path: Path,
    kind: str,
) -> None:
    release = metabolic_release(tmp_path)
    relation = release / "raw" / "reaction_ko.tsv"
    overlay = tmp_path / "overlay-reaction-ko.tsv"
    overlay.write_bytes(relation.read_bytes())
    relation.unlink()
    archive = _archive_release(release, tmp_path / f"overlay.{kind}", kind=kind)

    with pytest.raises(ValueError, match="missing roles.*reaction_ko"):
        KEGGDatabase.from_metabolic_files(archive)

    path = tmp_path / "archive-overlay.duckdb"
    KEGGDatabase.from_metabolic_files(archive, reaction_ko=overlay).write_duckdb(path)
    with duckdb.connect(str(path), read_only=True) as connection:
        source_paths = dict(
            connection.execute(
                "SELECT logical_name, display_path FROM _bioextract.source_file"
            ).fetchall()
        )
    assert source_paths["release_archive:0"] == str(archive)
    assert source_paths["reaction_ko:0"] == str(overlay)


def test_fully_replaced_archive_is_omitted_from_final_provenance(
    tmp_path: Path,
) -> None:
    release = metabolic_release(tmp_path)
    archive = _archive_release(release, tmp_path / "release.zip", kind="zip")
    raw = release / "raw"
    path = tmp_path / "fully-replaced.duckdb"
    KEGGDatabase.from_metabolic_files(
        archive,
        compound_list=raw / "compound" / "list.tsv",
        compound_entries=raw / "compound" / "entries",
        reaction_list=raw / "reaction" / "list.tsv",
        reaction_entries=raw / "reaction" / "entries",
        enzyme_list=raw / "enzyme" / "list.tsv",
        enzyme_entries=raw / "enzyme" / "entries",
        module_list=raw / "module" / "list.tsv",
        module_entries=raw / "module" / "entries",
        compound_pubchem=raw / "compound_pubchem.tsv",
        compound_reaction=raw / "compound_reaction.tsv",
        reaction_enzyme=raw / "reaction_enzyme.tsv",
        reaction_ko=raw / "reaction_ko.tsv",
        reaction_module=raw / "reaction_module.tsv",
        reaction_pathway=raw / "reaction_pathway.tsv",
        module_pathway=raw / "module_pathway.tsv",
    ).write_duckdb(path)
    with duckdb.connect(str(path), read_only=True) as connection:
        display_paths = {
            row[0]
            for row in connection.execute(
                "SELECT display_path FROM _bioextract.source_file"
            ).fetchall()
        }
    assert str(archive) not in display_paths
    assert str(raw / "reaction_ko.tsv") in display_paths


def test_archive_layout_keeps_optional_lists_when_required_roles_are_overlaid(
    tmp_path: Path,
) -> None:
    release = metabolic_release(tmp_path)
    for entries in release.glob("raw/*/entries"):
        shutil.rmtree(entries)
    for relation in release.glob("raw/*.tsv"):
        relation.unlink()
    archive = _archive_release(release, tmp_path / "lists-only.zip", kind="zip")
    files = metabolic_files(tmp_path)

    path = tmp_path / "fully-overlaid.duckdb"
    metabolic_source(files, source=archive).write_duckdb(path)
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute("SELECT count(*) FROM compound").fetchone() == (2,)


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
            "input_id": "3.6.1.3",
            "input_namespace": "ec",
            "entity_type": "reaction",
            "entity_id": "R00001",
            "match_type": "exact",
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
        "reaction_id"
    ].to_list() == ["R00001"]


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
    assert result.filter(result["module_id"] == "M00001")["is_complete"].item()

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
