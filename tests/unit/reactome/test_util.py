from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from bioextract.reactome.constant import COLS_MAPPING_RAW, SCHEMA_MAPPING_RAW
from bioextract.reactome.util import read_mapping_family_frame


def test_mapping_family_reader_preserves_literal_quotes_and_evidence_grain(
    tmp_path: Path,
) -> None:
    source = tmp_path / "mapping.tsv"
    row = (
        "P04637\tR-HSA-1\thttps://reactome.org/R-HSA-1\t"
        'A pathway with a "literal" quote\tTAS\tHomo sapiens'
    )
    source.write_text(
        "\n".join([row, row, row.replace("\tTAS\t", "\tIEA\t")]) + "\n",
        encoding="utf-8",
    )

    frame = read_mapping_family_frame(
        source,
        columns=COLS_MAPPING_RAW,
        schema=SCHEMA_MAPPING_RAW,
        context="mapping fixture",
    )

    assert frame.columns == COLS_MAPPING_RAW
    assert frame.height == 2
    assert frame.get_column("pathway_name").to_list() == [
        'A pathway with a "literal" quote',
        'A pathway with a "literal" quote',
    ]
    assert frame.get_column("evidence_code").to_list() == ["TAS", "IEA"]


def test_mapping_family_reader_supports_reaction_shaped_literal_fields(
    tmp_path: Path,
) -> None:
    source = tmp_path / "reaction-shaped.tsv"
    source.write_text(
        'P04637\tR-HSA-1\thttps://reactome.org/R-HSA-1\tEnzyme "A" event\tTAS\tHomo sapiens\n',
        encoding="utf-8",
    )
    columns = [
        "uniprot_id",
        "reactome_reaction_id",
        "reactome_url",
        "reaction_name",
        "evidence_code",
        "species",
    ]

    frame = read_mapping_family_frame(
        source,
        columns=columns,
        schema=dict.fromkeys(columns, pl.String),
        context="reaction-shaped fixture",
    )

    assert frame.get_column("reaction_name").to_list() == ['Enzyme "A" event']


@pytest.mark.parametrize(
    "content, message",
    [
        ("P04637\tR-HSA-1\turl\tname\tTAS\n", "exactly 6"),
        ("P04637\tR-HSA-1\turl\tname\tTAS\tHomo sapiens\textra\n", "exactly 6"),
        ("P04637\tR-HSA-1\turl\t\tTAS\tHomo sapiens\n", "empty required"),
    ],
)
def test_mapping_family_reader_rejects_ragged_or_empty_records(
    tmp_path: Path,
    content: str,
    message: str,
) -> None:
    source = tmp_path / "invalid.tsv"
    source.write_text(content, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        read_mapping_family_frame(
            source,
            columns=COLS_MAPPING_RAW,
            schema=SCHEMA_MAPPING_RAW,
            context="mapping fixture",
        )
