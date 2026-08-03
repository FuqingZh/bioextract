from __future__ import annotations

import polars as pl
import pytest

from bioextract._shared import create_group_input_frames

SCHEMA_GROUPS = {"GroupId": pl.String}
SCHEMA_GROUP_INPUT_IDS = {"GroupId": pl.String, "InputId": pl.String}


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

    assert df_groups.columns == ["GroupId"]
    assert df_group_membership.columns == ["GroupId", "InputId"]
    assert df_input_ids.columns == ["InputId"]
    assert df_groups.to_dicts() == [
        {"GroupId": "A"},
        {"GroupId": "B"},
        {"GroupId": "C"},
    ]
    assert df_group_membership.to_dicts() == [
        {"GroupId": "A", "InputId": "P04637"},
        {"GroupId": "A", "InputId": "TP53"},
        {"GroupId": "B", "InputId": "EGFR"},
        {"GroupId": "B", "InputId": "P04637"},
    ]
    assert df_input_ids.to_dicts() == [
        {"InputId": "EGFR"},
        {"InputId": "P04637"},
        {"InputId": "TP53"},
    ]


def test_create_group_input_frames_keeps_empty_group_membership_typed() -> None:
    group_input_frames = create_group_input_frames(
        {"A": []},
        schema_groups=SCHEMA_GROUPS,
        schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
    )

    assert group_input_frames.df_groups.to_dicts() == [{"GroupId": "A"}]
    assert group_input_frames.df_group_membership.schema == SCHEMA_GROUP_INPUT_IDS
    assert group_input_frames.df_group_membership.is_empty()
    assert group_input_frames.df_input_ids.schema == {"InputId": pl.String}
    assert group_input_frames.df_input_ids.is_empty()


def test_create_group_input_frames_rejects_duplicate_normalized_groups() -> None:
    with pytest.raises(ValueError, match="unique after normalization"):
        create_group_input_frames(
            {"A": ["TP53"], " A ": ["EGFR"]},
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
