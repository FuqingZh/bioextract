from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from bioextract.interpro import InterProDatabase


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
<release>
  <dbinfo dbname="INTERPRO" version="108.0"/>
</release>
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
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
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
    db = InterProDatabase.from_mapping_files(protein_to_interpro=files["protein2ipr"])

    row = (
        db.select_ids(
            ["P12345"],
        )
        .extract_mapping()
        .row(0, named=True)
    )

    assert row["InterProType"] is None
    assert row["MemberDb"] is None


def test_select_ids_streams_subset_and_reports_unmapped(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
    )

    selection = db.select_ids(
        ["sp|P12345|TEST_HUMAN", "MISSING"],
    )

    assert selection.extract_mapping().select(
        "InputId",
        "InputNamespace",
        "UniProtId",
        "InterProId",
        "MemberDb",
    ).to_dicts() == [
        {
            "InputId": "P12345",
            "InputNamespace": "uniprot",
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "MemberDb": "PFAM",
        },
        {
            "InputId": "P12345",
            "InputNamespace": "uniprot",
            "UniProtId": "P12345",
            "InterProId": "IPR000001",
            "MemberDb": "SMART",
        },
    ]
    assert selection.extract_unmatched_ids().to_dicts() == [{"InputId": "MISSING"}]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
    )

    grouped = db.select_groups(
        {"up": ["P12345"], "down": ["P12345", "MISSING"]},
    )

    df_mapping = grouped.extract_mapping()
    assert grouped.extract_mapping() is df_mapping
    assert df_mapping.columns[:3] == [
        "GroupId",
        "InputId",
        "InputNamespace",
    ]
    assert df_mapping.select("GroupId", "MemberDb").to_dicts() == [
        {"GroupId": "down", "MemberDb": "PFAM"},
        {"GroupId": "down", "MemberDb": "SMART"},
        {"GroupId": "up", "MemberDb": "PFAM"},
        {"GroupId": "up", "MemberDb": "SMART"},
    ]
    assert grouped.extract_unmatched_ids().to_dicts() == [
        {"GroupId": "down", "InputId": "MISSING"}
    ]


def test_write_parquet_writes_mapping_without_sidecar(tmp_path: Path) -> None:
    files = write_interpro_fixture(tmp_path)
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
    )

    path = tmp_path / "interpro.parquet"
    result = db.write_parquet(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    assert pl.read_parquet(path).height == 3


@pytest.mark.parametrize(
    ("mapping_row", "message"),
    [
        (
            "P12345\tIPR999999\tUnknown\tPF00051\t10\t80\n",
            "absent from XML entry metadata",
        ),
        (
            "P12345\tIPR000001\tKringle\tPF99999\t10\t80\n",
            "absent from XML member metadata",
        ),
    ],
)
def test_tidy_publication_rejects_mapping_relationships_absent_from_xml(
    tmp_path: Path,
    mapping_row: str,
    message: str,
) -> None:
    files = write_interpro_fixture(tmp_path)
    with gzip.open(files["protein2ipr"], "wt", encoding="utf-8") as handle:
        handle.write(mapping_row)
    database = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
    )

    with pytest.raises(ValueError, match=message):
        database.build_tidy()


@pytest.mark.parametrize(
    ("old_xml", "new_xml", "message"),
    [
        (
            'id="IPR000001" type="Domain"',
            'id="IPR000001" type=""',
            "one non-empty XML InterProType",
        ),
        (
            '<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>',
            (
                '<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>'
                '<db_xref db="OTHER" dbkey="PF00051" name="Kringle"/>'
            ),
            "one non-empty XML MemberDb",
        ),
        (
            "</member_list>\n</interpro>",
            (
                "</member_list>\n</interpro>"
                '<interpro id="IPR000001" type=""><name>Duplicate</name>'
                "</interpro>"
            ),
            "one non-empty XML InterProType",
        ),
        (
            '<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>',
            (
                '<db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>'
                '<db_xref db="" dbkey="PF00051" name="Kringle"/>'
            ),
            "one non-empty XML MemberDb",
        ),
    ],
)
def test_tidy_publication_requires_unique_nonempty_xml_enrichment(
    tmp_path: Path,
    old_xml: str,
    new_xml: str,
    message: str,
) -> None:
    files = write_interpro_fixture(tmp_path)
    with gzip.open(files["xml"], "rt", encoding="utf-8") as handle:
        xml = handle.read()
    with gzip.open(files["xml"], "wt", encoding="utf-8") as handle:
        handle.write(xml.replace(old_xml, new_xml, 1))
    database = InterProDatabase.from_mapping_files(
        protein_to_interpro=files["protein2ipr"],
        interpro_xml=files["xml"],
    )

    with pytest.raises(ValueError, match=message):
        database.build_tidy()
