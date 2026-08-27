from __future__ import annotations

import pytest

from bioextract._shared import normalize_uniprot_selection_id


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (" P04637 ", "P04637"),
        ("sp|P04637|P53_HUMAN", "P04637"),
        (" tr|A0A024RBG1|NUD4B_HUMAN ", "A0A024RBG1"),
        ("sp| P04637 | P53 HUMAN ", "P04637"),
        ("", ""),
        ("  ", ""),
    ],
)
def test_normalize_uniprot_selection_id(raw: str, expected: str) -> None:
    assert normalize_uniprot_selection_id(raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        "db|P04637|P53_HUMAN",
        "sp||P53_HUMAN",
        "sp|P04637|",
        "sp|P04637",
        "sp|P04637|P53_HUMAN|extra",
        "|P04637|P53_HUMAN",
        "sp|   |P53_HUMAN",
        "sp|P04637|   ",
    ],
)
def test_normalize_uniprot_selection_id_rejects_other_pipe_forms(raw: str) -> None:
    with pytest.raises(ValueError, match="Invalid UniProt pipe-form identifier"):
        normalize_uniprot_selection_id(raw)
