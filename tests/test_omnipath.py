from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from bioextract.omnipath import OmniPathDb, OmniPathResourceLimits


def _write_text_or_gzip(file_out: Path, content: str, *, should_gzip: bool) -> None:
    if should_gzip:
        with gzip.open(file_out, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return

    file_out.write_text(content, encoding="utf-8")


def _write_demo_omnipath_files(
    *,
    file_enzsub: Path | None = None,
    file_interactions: Path | None = None,
    should_gzip: bool = False,
) -> None:
    if file_enzsub is not None:
        _write_text_or_gzip(
            file_enzsub,
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

    if file_interactions is not None:
        _write_text_or_gzip(
            file_interactions,
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


def test_extract_enzsub_single_query_smoke(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(file_enzsub=file_enzsub)

    selection = OmniPathDb.from_files(file_enzsub=file_enzsub).select_ids(["P31749"])

    df_enzsub = selection.extract_enzsub()
    df_unmapped = selection.extract_unmapped_input_ids()

    assert df_enzsub.columns == ["SourceId", "TargetId", "TargetSite", "Modification"]
    assert df_enzsub.to_dicts() == [
        {
            "SourceId": "P31749",
            "TargetId": "BAD",
            "TargetSite": "S136",
            "Modification": "phosphorylation",
        },
        {
            "SourceId": "P31749",
            "TargetId": "FOXO3",
            "TargetSite": "T32",
            "Modification": "phosphorylation",
        },
    ]
    assert df_unmapped.to_dicts() == []


def test_extract_interactions_single_query_smoke(tmp_path: Path) -> None:
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(file_interactions=file_interactions)

    selection = (
        OmniPathDb.from_files(file_interactions=file_interactions)
        .select_ids(["ERBB2", "MISSING"])
        .with_interactions()
    )

    df_interactions = selection.extract_interactions()
    df_unmapped = selection.extract_unmapped_input_ids()

    assert df_interactions.columns == [
        "SourceId",
        "TargetId",
        "IsDirected",
        "IsStimulation",
        "IsInhibition",
    ]
    assert df_interactions.to_dicts() == [
        {
            "SourceId": "EGFR",
            "TargetId": "ERBB2",
            "IsDirected": False,
            "IsStimulation": True,
            "IsInhibition": False,
        }
    ]
    assert df_unmapped.to_dicts() == [{"InputId": "MISSING"}]


def test_group_selection_extracts_flat_outputs(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    selection = OmniPathDb.from_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    ).select_groups(
        {
            "G1": ["P31749", "ERBB2", "MISSING"],
            "G2": ["MAPK1", "MTOR"],
        }
    )

    df_enzsub = selection.extract_enzsub()
    df_interactions = selection.extract_interactions()
    df_unmapped = selection.extract_unmapped_input_ids()

    assert selection.is_grouped is True
    assert df_enzsub.columns == [
        "GroupId",
        "SourceId",
        "TargetId",
        "TargetSite",
        "Modification",
    ]
    assert df_interactions.columns == [
        "GroupId",
        "SourceId",
        "TargetId",
        "IsDirected",
        "IsStimulation",
        "IsInhibition",
    ]
    assert df_unmapped.columns == ["GroupId", "InputId"]
    assert df_enzsub.to_dicts() == [
        {
            "GroupId": "G1",
            "SourceId": "P31749",
            "TargetId": "BAD",
            "TargetSite": "S136",
            "Modification": "phosphorylation",
        },
        {
            "GroupId": "G1",
            "SourceId": "P31749",
            "TargetId": "FOXO3",
            "TargetSite": "T32",
            "Modification": "phosphorylation",
        },
        {
            "GroupId": "G2",
            "SourceId": "MAPK1",
            "TargetId": "ELK1",
            "TargetSite": "S383",
            "Modification": "phosphorylation",
        },
    ]
    assert df_interactions.to_dicts() == [
        {
            "GroupId": "G1",
            "SourceId": "EGFR",
            "TargetId": "ERBB2",
            "IsDirected": False,
            "IsStimulation": True,
            "IsInhibition": False,
        },
        {
            "GroupId": "G2",
            "SourceId": "AKT1",
            "TargetId": "MTOR",
            "IsDirected": True,
            "IsStimulation": True,
            "IsInhibition": False,
        },
        {
            "GroupId": "G2",
            "SourceId": "MTOR",
            "TargetId": "RPTOR",
            "IsDirected": True,
            "IsStimulation": True,
            "IsInhibition": False,
        },
    ]
    assert df_unmapped.to_dicts() == [{"GroupId": "G1", "InputId": "MISSING"}]


def test_extract_plain_and_gzip_inputs(tmp_path: Path) -> None:
    for should_gzip in (False, True):
        suffix = ".tsv.gz" if should_gzip else ".tsv"
        file_enzsub = tmp_path / f"enzsub{suffix}"
        file_interactions = tmp_path / f"interactions{suffix}"
        _write_demo_omnipath_files(
            file_enzsub=file_enzsub,
            file_interactions=file_interactions,
            should_gzip=should_gzip,
        )

        selection = OmniPathDb.from_files(
            file_enzsub=file_enzsub,
            file_interactions=file_interactions,
        ).select_ids(["P31749", "ERBB2"])

        assert selection.extract_enzsub().height == 2
        assert selection.extract_interactions().height == 1


def test_with_resources_constrains_extraction_and_reuses_cache(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    selection = OmniPathDb.from_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    ).select_ids(["P31749", "ERBB2"])

    df_enzsub = selection.extract_enzsub()
    selection_enzsub = selection.with_enzsub()

    assert selection_enzsub.extract_enzsub() is df_enzsub
    with pytest.raises(ValueError, match="not enabled"):
        selection_enzsub.extract_interactions()

    selection_interactions = selection.with_interactions()
    assert selection_interactions.extract_interactions().to_dicts() == [
        {
            "SourceId": "EGFR",
            "TargetId": "ERBB2",
            "IsDirected": False,
            "IsStimulation": True,
            "IsInhibition": False,
        }
    ]


def test_repeated_extraction_reuses_cached_frames(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    selection = OmniPathDb.from_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    ).select_groups({"G1": ["P31749"], "G2": ["ERBB2"]})

    df_enzsub_first = selection.extract_enzsub()
    df_enzsub_second = selection.extract_enzsub()
    df_inter_first = selection.extract_interactions()
    df_inter_second = selection.extract_interactions()
    df_unmapped_first = selection.extract_unmapped_input_ids()
    df_unmapped_second = selection.extract_unmapped_input_ids()

    assert df_enzsub_first is df_enzsub_second
    assert df_inter_first is df_inter_second
    assert df_unmapped_first is df_unmapped_second


def test_from_files_rejects_missing_or_oversized_files(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    with pytest.raises(ValueError, match="At least one OmniPath resource file"):
        OmniPathDb.from_files()

    with pytest.raises(FileNotFoundError, match="enzsub file not found"):
        OmniPathDb.from_files(file_enzsub=tmp_path / "missing.tsv")

    with pytest.raises(ValueError, match="enzsub file"):
        OmniPathDb.from_files(
            file_enzsub=file_enzsub,
            limits=OmniPathResourceLimits(file_enzsub_bytes_max=1),
        )

    with pytest.raises(ValueError, match="interactions file"):
        OmniPathDb.from_files(
            file_interactions=file_interactions,
            limits=OmniPathResourceLimits(file_interactions_bytes_max=1),
        )


def test_selection_limits_and_group_shape_are_enforced(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(file_enzsub=file_enzsub)

    db = OmniPathDb.from_files(
        file_enzsub=file_enzsub,
        limits=OmniPathResourceLimits(num_input_ids_max=1, num_groups_max=1),
    )

    with pytest.raises(ValueError, match="Normalized input ID count exceeds"):
        db.select_ids(["P31749", "MAPK1"])

    with pytest.raises(ValueError, match="non-empty"):
        db.select_groups({"  ": ["P31749"]})

    with pytest.raises(ValueError, match="unique after normalization"):
        db.select_groups({"A": ["P31749"], " A ": ["MAPK1"]})

    with pytest.raises(ValueError, match="Group count exceeds"):
        db.select_groups({"G1": ["P31749"], "G2": ["MAPK1"]})


def test_extract_rejects_missing_required_resource_file(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    selection_enzsub_only = OmniPathDb.from_files(file_enzsub=file_enzsub).select_ids(
        ["P31749"]
    )
    selection_inter_only = OmniPathDb.from_files(
        file_interactions=file_interactions
    ).select_ids(["ERBB2"])

    with pytest.raises(ValueError, match="without interactions file"):
        selection_enzsub_only.with_interactions().extract_interactions()
    with pytest.raises(ValueError, match="without enzsub file"):
        selection_inter_only.with_enzsub().extract_enzsub()


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

    selection = OmniPathDb.from_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    ).select_ids(["P31749", "AKT1"])

    with pytest.raises(ValueError, match="enzsub file"):
        selection.extract_enzsub()
    with pytest.raises(ValueError, match="interactions file"):
        selection.extract_interactions()


def test_extract_unmapped_with_single_selected_resource(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    file_interactions = tmp_path / "interactions.tsv"
    _write_demo_omnipath_files(
        file_enzsub=file_enzsub,
        file_interactions=file_interactions,
    )

    selection = (
        OmniPathDb.from_files(
            file_enzsub=file_enzsub,
            file_interactions=file_interactions,
        )
        .select_ids(["P31749", "ERBB2"])
        .with_enzsub()
    )

    assert selection.extract_unmapped_input_ids().to_dicts() == [
        {"InputId": "ERBB2"}
    ]


def test_selection_exposes_group_mode(tmp_path: Path) -> None:
    file_enzsub = tmp_path / "enzsub.tsv"
    _write_demo_omnipath_files(file_enzsub=file_enzsub)
    db = OmniPathDb.from_files(file_enzsub=file_enzsub)

    assert db.select_ids(["P31749"]).is_grouped is False
    assert db.select_groups({"G1": ["P31749"]}).is_grouped is True
