"""Official-format fixtures exercise complex membership and integrity boundaries."""

import csv
from collections.abc import Iterator
from pathlib import Path

import polars as pl
import pytest

from bioextract import ComplexPortalDatabase
from bioextract.complexportal.complexportal import COLUMNS


def _file(
    tmp_path: Path, members: str = "P12345(2)|[Q12345,Q12346](0)|P12346-PRO_00001(1)"
) -> Path:
    path = tmp_path / "arbitrary-filename.tsv"
    row = dict.fromkeys(COLUMNS, "-")
    row.update(
        {
            "#Complex ac": "CPX-17",
            "Recommended name": "Known complex",
            "Taxonomy identifier": "9606",
            "Expanded participant list": members,
            "Evidence Code": "ECO:0000353(physical interaction)",
            "Identifiers (and stoichiometry) of molecules in complex": "P12345(2)|CPX-1(1)",
            "Description": "Fixture",
            "Source": "IntAct",
        }
    )
    with path.open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=COLUMNS, delimiter="\t", quoting=csv.QUOTE_NONE
        )
        writer.writeheader()
        writer.writerow(row)
    return path


def test_exact_source_columns_alternatives_and_unknown_stoichiometry(
    tmp_path: Path,
) -> None:
    path = _file(tmp_path)
    db = ComplexPortalDatabase.from_files(complextab=path, release_version="2026-01-09")
    assert tuple(db.complexes().collect_schema().names()) == COLUMNS
    members = db.expanded_members().collect()
    assert members.height == 4
    assert members.filter(pl.col("participant_id") == "P12345")["stoichiometry"][0] == 2
    alternatives = members.filter(pl.col("is_alternative"))
    assert alternatives.height == 2 and alternatives["member_group_id"].n_unique() == 1
    assert alternatives["stoichiometry"].null_count() == 2
    assert "P12346-PRO_00001" in members["participant_id"]
    empty = db.with_taxon("39947").expanded_members().collect()
    assert empty.is_empty() and empty.schema == members.schema
    assert db.with_taxon("9606").expanded_members().collect().equals(members)


@pytest.mark.parametrize(
    "token",
    [
        "P12345",
        "P12345(-1)",
        "P12345(1.5)",
        "CPX-1(1)",
        "[P12345,](1)",
        "P12345(1)|P12345(1)",
    ],
)
def test_malformed_or_duplicated_protein_tokens_fail(
    tmp_path: Path, token: str
) -> None:
    with pytest.raises(ValueError, match="token"):
        ComplexPortalDatabase.from_files(
            complextab=_file(tmp_path, token), release_version="2026-01-09"
        )


def test_invalid_header_missing_release_and_source_mutation_fail(
    tmp_path: Path,
) -> None:
    path = _file(tmp_path)
    with pytest.raises(ValueError, match="explicit"):
        ComplexPortalDatabase.from_files(complextab=path, release_version="")
    db = ComplexPortalDatabase.from_files(complextab=path, release_version="2026-01-09")
    path.write_text(path.read_text().replace("CPX-17", "CPX-18"))
    with pytest.raises(ValueError, match="changed"):
        db.expanded_members()
    path.write_text(path.read_text().replace("#Complex ac", "complex_id"))
    with pytest.raises(ValueError, match="header"):
        ComplexPortalDatabase.from_files(complextab=path, release_version="2026-01-09")


def test_unresolved_intact_participant_is_preserved_without_uniprot_inference(
    tmp_path: Path,
) -> None:
    db = ComplexPortalDatabase.from_files(
        complextab=_file(tmp_path, "EBI-25775462(2)"), release_version="2026-01-09"
    )
    row = db.expanded_members().collect().row(0, named=True)
    assert row["participant_id"] == "EBI-25775462"
    assert row["identifier_namespace"] == "intact"


@pytest.mark.parametrize("accession", ["ABCDEF", "P123456", "A12345", "A0A023GPIX"])
def test_invalid_uniprot_accession_is_rejected(tmp_path: Path, accession: str) -> None:
    with pytest.raises(ValueError, match="token"):
        ComplexPortalDatabase.from_files(
            complextab=_file(tmp_path, f"{accession}(1)"), release_version="test"
        )


@pytest.mark.parametrize("accession", ["A2BC19", "P12345", "A0A023GPI8", "Q12345-2"])
def test_official_uniprot_accession_forms_are_retained(
    tmp_path: Path, accession: str
) -> None:
    db = ComplexPortalDatabase.from_files(
        complextab=_file(tmp_path, f"{accession}(1)"), release_version="test"
    )
    assert db.expanded_members().collect()["participant_id"].to_list() == [accession]


@pytest.mark.parametrize(
    "members",
    ["P12345(1)|P12345(2)", "[P12345,P12345](1)", "P12345(1)|[P12345,Q12345](1)"],
)
def test_duplicate_participant_identity_is_rejected(
    tmp_path: Path, members: str
) -> None:
    with pytest.raises(ValueError, match="duplicate"):
        ComplexPortalDatabase.from_files(
            complextab=_file(tmp_path, members), release_version="test"
        )


def test_stoichiometry_range_is_validated_before_returning_handle(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="stoichiometry"):
        ComplexPortalDatabase.from_files(
            complextab=_file(tmp_path, "P12345(9223372036854775808)"),
            release_version="test",
        )
    db = ComplexPortalDatabase.from_files(
        complextab=_file(tmp_path, "P12345(9223372036854775807)"),
        release_version="test",
    )
    assert db.expanded_members().collect()["stoichiometry"][0] == 9223372036854775807


@pytest.mark.parametrize("relation", ["complexes", "expanded_members"])
@pytest.mark.parametrize("replay", [False, True])
def test_lazy_execution_rejects_changed_source(
    tmp_path: Path, relation: str, replay: bool
) -> None:
    path = _file(tmp_path)
    db = ComplexPortalDatabase.from_files(complextab=path, release_version="test")
    plan = getattr(db.with_taxon("9606"), relation)().head(2)
    if replay:
        assert plan.collect().height > 0
    path.write_text(path.read_text().replace("CPX-17", "CPX-18"))
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="ComplexTab changed"
    ):
        plan.collect()


def test_replay_preserves_filter_projection_and_limit(tmp_path: Path) -> None:
    path = _file(tmp_path)
    text = path.read_text()
    row = text.splitlines()[-1].replace("CPX-17", "CPX-18").replace("9606", "10090")
    path.write_text(text + row + "\n")
    db = ComplexPortalDatabase.from_files(complextab=path, release_version="test")
    plan = (
        db.complexes()
        .filter(pl.col("Taxonomy identifier") == "10090")
        .select("#Complex ac")
        .head(1)
    )
    assert plan.collect().to_dicts() == [{"#Complex ac": "CPX-18"}]
    assert plan.collect().to_dicts() == [{"#Complex ac": "CPX-18"}]


def test_source_change_during_native_scan_cannot_emit_new_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _file(tmp_path)
    db = ComplexPortalDatabase.from_files(complextab=path, release_version="test")
    plan = db.complexes()
    original_batches = pl.LazyFrame.collect_batches

    def changed_batches(
        frame: pl.LazyFrame, *, chunk_size: int
    ) -> Iterator[pl.DataFrame]:
        for batch in original_batches(frame, chunk_size=chunk_size):
            path.write_text(path.read_text().replace("CPX-17", "CPX-18"))
            yield batch

    monkeypatch.setattr(pl.LazyFrame, "collect_batches", changed_batches)
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="ComplexTab changed"
    ):
        plan.collect()
