from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

import bioextract.omnipath.util as omnipath_util
from bioextract.omnipath import OmniPathDatabase


def _write_text_or_gzip(path: Path, content: str, *, should_gzip: bool) -> None:
    if should_gzip:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return

    path.write_text(content, encoding="utf-8")


def _write_demo_omnipath_files(
    *,
    enzsub: Path | None = None,
    interactions: Path | None = None,
    should_gzip: bool = False,
) -> None:
    if enzsub is not None:
        _write_text_or_gzip(
            enzsub,
            "\n".join(
                [
                    "enzyme\tsubstrate\tresidue_type\tresidue_offset\tmodification",
                    "P31749\tBAD\tS\t136\tphosphorylation",
                    "P31749\tFOXO3\tT\t32\tphosphorylation",
                    "MAPK1\tELK1\tS\t383\tphosphorylation",
                ]
            )
            + "\n",
            should_gzip=should_gzip,
        )

    if interactions is not None:
        _write_text_or_gzip(
            interactions,
            "\n".join(
                [
                    "source\ttarget\tis_directed\tis_stimulation\tis_inhibition",
                    "AKT1\tMTOR\tTrue\tTrue\tFalse",
                    "MTOR\tRPTOR\tTrue\tTrue\tFalse",
                    "EGFR\tERBB2\tFalse\tTrue\tFalse",
                ]
            )
            + "\n",
            should_gzip=should_gzip,
        )


def test_extract_enzsub_minimal_round_trip(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(enzsub=file_enzsub)

    selection = OmniPathDatabase.from_files(enzsub=file_enzsub).select_ids(["P31749"])

    df_enzsub = selection.enzsub().collect()
    df_unmapped = selection.unmatched_ids().collect()

    assert df_enzsub.columns == [
        "enzyme",
        "substrate",
        "residue_type",
        "residue_offset",
        "modification",
        "target_site",
    ]
    assert df_enzsub.to_dicts() == [
        {
            "enzyme": "P31749",
            "substrate": "BAD",
            "residue_type": "S",
            "residue_offset": "136",
            "modification": "phosphorylation",
            "target_site": "S136",
        },
        {
            "enzyme": "P31749",
            "substrate": "FOXO3",
            "residue_type": "T",
            "residue_offset": "32",
            "modification": "phosphorylation",
            "target_site": "T32",
        },
    ]
    assert df_unmapped.to_dicts() == []


def test_extract_interactions_minimal_round_trip(tmp_path: Path) -> None:
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(interactions=file_interactions)

    selection = (
        OmniPathDatabase.from_files(interactions=file_interactions)
        .select_ids(["ERBB2", "MISSING"])
        .with_interactions()
    )

    df_interactions = selection.interactions().collect()
    df_unmapped = selection.unmatched_ids().collect()

    assert df_interactions.columns == [
        "source",
        "target",
        "is_directed",
        "is_stimulation",
        "is_inhibition",
    ]
    assert df_interactions.to_dicts() == [
        {
            "source": "EGFR",
            "target": "ERBB2",
            "is_directed": False,
            "is_stimulation": True,
            "is_inhibition": False,
        }
    ]
    assert df_unmapped.to_dicts() == [{"input_id": "MISSING"}]


def test_group_selection_extracts_flat_outputs(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    selection = OmniPathDatabase.from_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    ).select_groups(
        {
            "G1": ["P31749", "ERBB2", "MISSING"],
            "G2": ["MAPK1", "MTOR"],
        }
    )

    df_enzsub = selection.enzsub().collect()
    df_interactions = selection.interactions().collect()
    df_unmapped = selection.unmatched_ids().collect()

    assert selection.is_grouped is True
    assert df_enzsub.columns == [
        "group_id",
        "enzyme",
        "substrate",
        "residue_type",
        "residue_offset",
        "modification",
        "target_site",
    ]
    assert df_interactions.columns == [
        "group_id",
        "source",
        "target",
        "is_directed",
        "is_stimulation",
        "is_inhibition",
    ]
    assert df_unmapped.columns == ["group_id", "input_id"]
    assert df_enzsub.to_dicts() == [
        {
            "group_id": "G1",
            "enzyme": "P31749",
            "substrate": "BAD",
            "residue_type": "S",
            "residue_offset": "136",
            "modification": "phosphorylation",
            "target_site": "S136",
        },
        {
            "group_id": "G1",
            "enzyme": "P31749",
            "substrate": "FOXO3",
            "residue_type": "T",
            "residue_offset": "32",
            "modification": "phosphorylation",
            "target_site": "T32",
        },
        {
            "group_id": "G2",
            "enzyme": "MAPK1",
            "substrate": "ELK1",
            "residue_type": "S",
            "residue_offset": "383",
            "modification": "phosphorylation",
            "target_site": "S383",
        },
    ]
    assert df_interactions.to_dicts() == [
        {
            "group_id": "G1",
            "source": "EGFR",
            "target": "ERBB2",
            "is_directed": False,
            "is_stimulation": True,
            "is_inhibition": False,
        },
        {
            "group_id": "G2",
            "source": "AKT1",
            "target": "MTOR",
            "is_directed": True,
            "is_stimulation": True,
            "is_inhibition": False,
        },
        {
            "group_id": "G2",
            "source": "MTOR",
            "target": "RPTOR",
            "is_directed": True,
            "is_stimulation": True,
            "is_inhibition": False,
        },
    ]
    assert df_unmapped.to_dicts() == [{"group_id": "G1", "input_id": "MISSING"}]


def test_extract_plain_and_gzip_inputs(tmp_path: Path) -> None:
    for should_gzip in (False, True):
        suffix = ".tsv.gz" if should_gzip else ".tsv"
        file_enzsub = tmp_path / f"enzsub{suffix}"
        file_interactions = tmp_path / f"interactions{suffix}"
        _write_demo_omnipath_files(
            enzsub=file_enzsub,
            interactions=file_interactions,
            should_gzip=should_gzip,
        )

        selection = OmniPathDatabase.from_files(
            enzsub=file_enzsub,
            interactions=file_interactions,
        ).select_ids(["P31749", "ERBB2"])

        assert selection.enzsub().collect().height == 2
        assert selection.interactions().collect().height == 1


def test_source_scans_preserve_official_columns(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(enzsub=file_enzsub, interactions=file_interactions)

    db = OmniPathDatabase.from_files(enzsub=file_enzsub, interactions=file_interactions)
    assert db.scan_enzsub().collect_schema().names() == [
        "enzyme",
        "substrate",
        "residue_type",
        "residue_offset",
        "modification",
    ]
    assert db.scan_interactions().collect_schema().names() == [
        "source",
        "target",
        "is_directed",
        "is_stimulation",
        "is_inhibition",
    ]


def test_with_resources_constrains_extraction_and_reuses_cache(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    selection = OmniPathDatabase.from_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    ).select_ids(["P31749", "ERBB2"])

    df_enzsub = selection.enzsub().collect()
    selection_enzsub = selection.with_enzsub()

    assert selection_enzsub.enzsub().collect().equals(df_enzsub)
    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="not enabled"):
        selection_enzsub.interactions().collect()

    selection_interactions = selection.with_interactions()
    assert selection_interactions.interactions().collect().to_dicts() == [
        {
            "source": "EGFR",
            "target": "ERBB2",
            "is_directed": False,
            "is_stimulation": True,
            "is_inhibition": False,
        }
    ]


def test_repeated_extraction_reuses_cached_frames(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    selection = OmniPathDatabase.from_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    ).select_groups({"G1": ["P31749"], "G2": ["ERBB2"]})

    df_enzsub_first = selection.enzsub().collect()
    df_enzsub_second = selection.enzsub().collect()
    df_inter_first = selection.interactions().collect()
    df_inter_second = selection.interactions().collect()
    df_unmapped_first = selection.unmatched_ids().collect()
    df_unmapped_second = selection.unmatched_ids().collect()

    assert df_enzsub_first.equals(df_enzsub_second)
    assert df_inter_first.equals(df_inter_second)
    assert df_unmapped_first.equals(df_unmapped_second)


def test_from_files_rejects_missing_files(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="At least one OmniPath resource file"):
        OmniPathDatabase.from_files()

    with pytest.raises(FileNotFoundError, match="enzsub file not found"):
        OmniPathDatabase.from_files(enzsub=tmp_path / "missing.tsv")


def test_group_shape_is_enforced(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(enzsub=file_enzsub)

    db = OmniPathDatabase.from_files(enzsub=file_enzsub)

    with pytest.raises(ValueError, match="non-empty"):
        db.select_groups({"  ": ["P31749"]})

    with pytest.raises(ValueError, match="unique after normalization"):
        db.select_groups({"A": ["P31749"], " A ": ["MAPK1"]})


def test_extract_rejects_missing_required_resource_file(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    selection_enzsub_only = OmniPathDatabase.from_files(enzsub=file_enzsub).select_ids(
        ["P31749"]
    )
    selection_inter_only = OmniPathDatabase.from_files(
        interactions=file_interactions
    ).select_ids(["ERBB2"])

    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="without interactions file"
    ):
        selection_enzsub_only.with_interactions().interactions().collect()
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="without enzsub file"
    ):
        selection_inter_only.with_enzsub().enzsub().collect()


def test_extract_rejects_missing_required_columns(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_text_or_gzip(
        file_enzsub,
        "enzyme\tsubstrate\tmodification\nP31749\tBAD\tphosphorylation\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_interactions,
        "source\ttarget\tis_directed\nAKT1\tMTOR\tTrue\n",
        should_gzip=False,
    )

    with pytest.raises(ValueError, match="enzsub file"):
        OmniPathDatabase.from_files(
            enzsub=file_enzsub,
            interactions=file_interactions,
        )


def test_extract_unmapped_with_single_selected_resource(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    selection = (
        OmniPathDatabase.from_files(
            enzsub=file_enzsub,
            interactions=file_interactions,
        )
        .select_ids(["P31749", "ERBB2"])
        .with_enzsub()
    )

    assert selection.unmatched_ids().collect().to_dicts() == [{"input_id": "ERBB2"}]


def test_selection_exposes_group_mode(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(enzsub=file_enzsub)
    db = OmniPathDatabase.from_files(enzsub=file_enzsub)

    assert db.select_ids(["P31749"]).is_grouped is False
    assert db.select_groups({"G1": ["P31749"]}).is_grouped is True


def test_group_selection_matches_unique_ids_then_expands_membership(
    tmp_path: Path,
) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(enzsub=file_enzsub)

    selection = OmniPathDatabase.from_files(enzsub=file_enzsub).select_groups(
        {
            "G1": ["sp|P31749|AKT1_HUMAN", "MISSING"],
            "G2": ["P31749", "MISSING"],
        }
    )

    assert selection._df_input_ids.to_dicts() == [  # pyright: ignore[reportPrivateUsage]
        {"input_id": "MISSING"},
        {"input_id": "P31749"},
    ]
    group_membership = selection._df_group_membership  # pyright: ignore[reportPrivateUsage]
    assert group_membership is not None
    assert group_membership.to_dicts() == [
        {"group_id": "G1", "input_id": "MISSING"},
        {"group_id": "G1", "input_id": "P31749"},
        {"group_id": "G2", "input_id": "MISSING"},
        {"group_id": "G2", "input_id": "P31749"},
    ]
    assert (
        selection.enzsub()
        .collect()
        .group_by("group_id")
        .len()
        .sort("group_id")
        .to_dicts()
    ) == [
        {"group_id": "G1", "len": 2},
        {"group_id": "G2", "len": 2},
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "G1", "input_id": "MISSING"},
        {"group_id": "G2", "input_id": "MISSING"},
    ]


def test_group_selection_builds_one_scan_per_cached_resource(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    )

    original_scan_enzsub = omnipath_util.scan_enzsub
    original_scan_interactions = omnipath_util.scan_interactions
    scan_counts = {"enzsub": 0, "interactions": 0}

    def counted_scan_enzsub(file_enzsub: Path) -> pl.LazyFrame:
        scan_counts["enzsub"] += 1
        return original_scan_enzsub(file_enzsub)

    def counted_scan_interactions(file_interactions: Path) -> pl.LazyFrame:
        scan_counts["interactions"] += 1
        return original_scan_interactions(file_interactions)

    monkeypatch.setattr(omnipath_util, "scan_enzsub", counted_scan_enzsub)
    monkeypatch.setattr(
        omnipath_util,
        "scan_interactions",
        counted_scan_interactions,
    )

    selection = OmniPathDatabase.from_files(
        enzsub=file_enzsub,
        interactions=file_interactions,
    ).select_groups(
        {
            "G1": ["sp|P31749|AKT1_HUMAN", "ERBB2"],
            "G2": ["P31749", "ERBB2"],
        }
    )
    assert selection.enzsub().collect().equals(selection.enzsub().collect())
    assert selection.interactions().collect().equals(selection.interactions().collect())
    assert (
        selection.unmatched_ids().collect().equals(selection.unmatched_ids().collect())
    )
    assert all(value >= 1 for value in scan_counts.values())
