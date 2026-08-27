from __future__ import annotations

import polars as pl
import pytest

from bioextract._shared import create_group_input_frames, create_input_id_frame

SCHEMA_GROUPS = {"group_id": pl.String}
SCHEMA_GROUP_INPUT_IDS = {"group_id": pl.String, "input_id": pl.String}
SCHEMA_UNMAPPED = {"input_id": pl.String}


def test_create_input_id_frame_is_trim_only() -> None:
    frame = create_input_id_frame(
        [" sp|P04637|P53_HUMAN ", "P04637", ""],
        schema_unmapped=SCHEMA_UNMAPPED,
    )

    assert frame.to_dicts() == [
        {"input_id": "P04637"},
        {"input_id": "sp|P04637|P53_HUMAN"},
    ]


def test_create_group_input_frames_preserves_group_contract() -> None:
    group_input_frames = create_group_input_frames(
        {
            " B ": [" EGFR ", " ", "EGFR", "sp|P04637|P53_HUMAN"],
            "A": ["sp|P04637|P53_HUMAN", "TP53", "P04637"],
            "C": [],
        },
        schema_groups=SCHEMA_GROUPS,
        schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
    )
    df_groups = group_input_frames.df_groups
    df_group_membership = group_input_frames.df_group_membership
    df_input_ids = group_input_frames.df_input_ids

    assert df_groups.columns == ["group_id"]
    assert df_group_membership.columns == ["group_id", "input_id"]
    assert df_input_ids.columns == ["input_id"]
    assert df_groups.to_dicts() == [
        {"group_id": "A"},
        {"group_id": "B"},
        {"group_id": "C"},
    ]
    assert df_group_membership.to_dicts() == [
        {"group_id": "A", "input_id": "P04637"},
        {"group_id": "A", "input_id": "TP53"},
        {"group_id": "A", "input_id": "sp|P04637|P53_HUMAN"},
        {"group_id": "B", "input_id": "EGFR"},
        {"group_id": "B", "input_id": "sp|P04637|P53_HUMAN"},
    ]
    assert df_input_ids.to_dicts() == [
        {"input_id": "EGFR"},
        {"input_id": "P04637"},
        {"input_id": "TP53"},
        {"input_id": "sp|P04637|P53_HUMAN"},
    ]


def test_create_group_input_frames_keeps_empty_group_membership_typed() -> None:
    group_input_frames = create_group_input_frames(
        {"A": []},
        schema_groups=SCHEMA_GROUPS,
        schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
    )

    assert group_input_frames.df_groups.to_dicts() == [{"group_id": "A"}]
    assert group_input_frames.df_group_membership.schema == SCHEMA_GROUP_INPUT_IDS
    assert group_input_frames.df_group_membership.is_empty()
    assert group_input_frames.df_input_ids.schema == {"input_id": pl.String}
    assert group_input_frames.df_input_ids.is_empty()


def test_create_group_input_frames_rejects_duplicate_normalized_groups() -> None:
    with pytest.raises(ValueError, match="unique after normalization"):
        create_group_input_frames(
            {"A": ["TP53"], " A ": ["EGFR"]},
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
