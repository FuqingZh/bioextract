from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from bioextract.interpro import InterProDb


def write_interpro_fixture(tmp_path: Path) -> dict[str, Path]:
    file_protein2ipr = tmp_path / "protein2ipr.dat.gz"
    with gzip.open(file_protein2ipr, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "P12345\tIPR000001\tKringle\tPF00051\t10\t80",
                    "P12345\tIPR000001\tKringle\tSM00130\t12\t76",
                    "Q9Y243\tIPR000002\tProtein kinase domain\tPF00069\t5\t220",
                ]
            )
            + "\n"
        )

    file_xml = tmp_path / "interpro.xml.gz"
    with gzip.open(file_xml, "wt", encoding="utf-8") as handle:
        handle.write(
            """<?xml version="1.0" encoding="UTF-8"?>
<interprodb>
<interpro id="IPR000001" type="Domain">
  <name>Kringle</name>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>
    <db_xref db="SMART" dbkey="SM00130" name="KR"/>
  </member_list>
</interpro>
<interpro id="IPR000002" type="Homologous_superfamily">
  <name>Protein kinase domain</name>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00069" name="Pkinase"/>
  </member_list>
</interpro>
</interprodb>
"""
        )
    return {"protein2ipr": file_protein2ipr, "xml": file_xml}


def test_extract_mapping_uses_xml_metadata_when_available(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=files["protein2ipr"],
        file_interpro_xml=files["xml"],
    )

    assert db.extract_mapping().to_dicts() == [
        {
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "InterProName": "Kringle",
            "InterProType": "Domain",
            "MemberDb": "PFAM",
            "MemberDbId": "PF00051",
            "Start": 10,
            "End": 80,
        },
        {
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "InterProName": "Kringle",
            "InterProType": "Domain",
            "MemberDb": "SMART",
            "MemberDbId": "SM00130",
            "Start": 12,
            "End": 76,
        },
        {
            "UniProtId": "Q9Y243",
            "InterProId": "IPR000002",
            "InterProName": "Protein kinase domain",
            "InterProType": "Homologous_superfamily",
            "MemberDb": "PFAM",
            "MemberDbId": "PF00069",
            "Start": 5,
            "End": 220,
        },
    ]


def test_xml_metadata_is_optional(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(file_protein2ipr=files["protein2ipr"])

    row = (
        db.select_ids(
            ["P12345"],
            kind_input_id="uniprot",
        )
        .extract_mapping()
        .row(0, named=True)
    )

    assert row["InterProType"] is None
    assert row["MemberDb"] is None


def test_select_ids_streams_subset_and_reports_unmapped(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=files["protein2ipr"],
        file_interpro_xml=files["xml"],
    )

    selection = db.select_ids(
        ["sp|P12345|TEST_HUMAN", "MISSING"],
        kind_input_id="uniprot",
    )

    assert selection.extract_mapping().select(
        "InputId",
        "KindInputId",
        "UniProtId",
        "InterProId",
        "MemberDb",
    ).to_dicts() == [
        {
            "InputId": "P12345",
            "KindInputId": "uniprot",
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "MemberDb": "PFAM",
        },
        {
            "InputId": "P12345",
            "KindInputId": "uniprot",
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "MemberDb": "SMART",
        },
    ]
    assert selection.extract_unmapped_input_ids().to_dicts() == [{"InputId": "MISSING"}]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=files["protein2ipr"],
        file_interpro_xml=files["xml"],
    )

    grouped = db.select_groups(
        {"up": ["P12345"], "down": ["MISSING"]},
        kind_input_id="uniprot",
    )

    assert grouped.extract_mapping().columns[:3] == [
        "GroupId",
        "InputId",
        "KindInputId",
    ]
    assert grouped.extract_mapping().select("GroupId", "MemberDb").to_dicts() == [
        {"GroupId": "up", "MemberDb": "PFAM"},
        {"GroupId": "up", "MemberDb": "SMART"},
    ]
    assert grouped.extract_unmapped_input_ids().to_dicts() == [
        {"GroupId": "down", "InputId": "MISSING"}
    ]


def test_build_tidy_writes_mapping_parquet_and_manifest(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=files["protein2ipr"],
        file_interpro_xml=files["xml"],
    )

    report = db.write_tidy(tmp_path / "out", should_write_manifest=True)

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "interpro-mapping-v0.1"
    assert [asset.path for asset in report.assets] == ["mapping.parquet"]
    assert pl.read_parquet(tmp_path / "out" / "mapping.parquet").height == 3


def test_validates_kind_input_id(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDb.from_mapping_files(file_protein2ipr=files["protein2ipr"])

    with pytest.raises(ValueError, match="kind_input_id"):
        db.select_ids(["P12345"], kind_input_id="geneid")  # type: ignore[arg-type]
