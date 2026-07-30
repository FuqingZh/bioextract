from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import polars as pl
import pytest
from polars.testing import assert_frame_equal

import bioextract.interpro as interpro
from bioextract.interpro import InterProDatabase


def write_interpro_snapshot(
    dir_parent: Path,
    *,
    version: str = "108.0",
    protein2ipr_rows: list[str] | None = None,
    xml_entries: str | None = None,
) -> tuple[Path, Path]:
    dir_raw = dir_parent / version / "raw"
    dir_raw.mkdir(parents=True)
    file_protein2ipr = dir_raw / "protein2ipr.dat.gz"
    with gzip.open(file_protein2ipr, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                protein2ipr_rows
                or [
                    "P12345\tIPR000001\tInterPro Kringle domain\tPF00051\t10\t80",
                    "P12345\tIPR000001\tInterPro Kringle domain\tPF00051\t100\t180",
                    "P12345\tIPR000002\tProtein kinase domain\tPF00069\t5\t220",
                    "Q9Y243\tIPR000002\tProtein kinase domain\tPF00069\t8\t210",
                    "P12345\tIPR000001\tInterPro Kringle domain\tSM00130\t12\t76",
                ]
            )
            + "\n"
        )

    file_xml = dir_raw / "interpro.xml.gz"
    entries = (
        xml_entries
        or """
<interpro id="IPR000001" type="Domain">
  <name>InterPro Kringle domain</name>
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
"""
    )
    with gzip.open(file_xml, "wt", encoding="utf-8") as handle:
        handle.write(
            f"""<?xml version="1.0" encoding="UTF-8"?>
<interprodb>
<release>
  <dbinfo version="{version}" dbname="INTERPRO"/>
</release>
{entries}
</interprodb>
"""
        )
    return file_protein2ipr, file_xml


def test_write_pfam_duckdb_emits_compact_relations(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    path = tmp_path / "interpro_pfam.duckdb"

    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=file_protein2ipr,
        interpro_xml=file_xml,
    )
    result = db.write_duckdb(path)

    assert result.tables == ("protein_term", "term", "term_xref")
    with duckdb.connect(str(path), read_only=True) as connection:
        df_protein_term = pl.read_database(  # pyright: ignore[reportUnknownMemberType]  # Polars-DuckDB boundary
            "SELECT * FROM protein_term", connection
        ).sort("uniprot_id", "pfam_id")
        df_term = pl.read_database(  # pyright: ignore[reportUnknownMemberType]  # Polars-DuckDB boundary
            "SELECT * FROM term", connection
        )
        df_term_xref = pl.read_database(  # pyright: ignore[reportUnknownMemberType]  # Polars-DuckDB boundary
            "SELECT * FROM term_xref", connection
        )
    assert df_protein_term.schema == {
        "uniprot_id": pl.String,
        "pfam_id": pl.String,
    }
    assert df_protein_term.to_dicts() == [
        {"uniprot_id": "P12345", "pfam_id": "PF00051"},
        {"uniprot_id": "P12345", "pfam_id": "PF00069"},
        {"uniprot_id": "Q9Y243", "pfam_id": "PF00069"},
    ]

    assert df_term.to_dicts() == [
        {"pfam_id": "PF00051", "pfam_name": "Kringle"},
        {"pfam_id": "PF00069", "pfam_name": "Pkinase"},
    ]
    assert "InterPro Kringle domain" not in df_term.get_column("pfam_name").to_list()

    assert df_term_xref.to_dicts() == [
        {
            "pfam_id": "PF00051",
            "interpro_id": "IPR000001",
            "interpro_name": "InterPro Kringle domain",
            "interpro_type": "Domain",
        },
        {
            "pfam_id": "PF00069",
            "interpro_id": "IPR000002",
            "interpro_name": "Protein kinase domain",
            "interpro_type": "Homologous_superfamily",
        },
    ]


def test_build_pfam_tidy_keeps_lazy_frames(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)

    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=file_protein2ipr,
        interpro_xml=file_xml,
    )
    dataset = db.build_tidy(config="pfam")
    result = dataset.write_duckdb(tmp_path / "interpro_pfam.duckdb")

    assert dataset.resource_schema_version == "interpro-pfam-v0.1"
    assert result.tables == ("protein_term", "term", "term_xref")
    assert all(isinstance(frame, pl.LazyFrame) for frame in dataset.frames.values())
    assert dataset.frames["protein_term"].collect_schema() == pl.Schema(
        {"UniProtId": pl.String, "PfamId": pl.String}
    )
    assert dataset.frames["term"].collect_schema() == pl.Schema(
        {"PfamId": pl.String, "PfamName": pl.String}
    )
    assert dataset.frames["term_xref"].collect_schema() == pl.Schema(
        {
            "PfamId": pl.String,
            "InterProId": pl.String,
            "InterProName": pl.String,
            "InterProType": pl.String,
        }
    )
    assert "maintain_order: true" not in dataset.frames["protein_term"].explain(
        engine="streaming"
    )


@pytest.mark.parametrize(
    ("protein2ipr_rows", "xml_entries", "error"),
    [
        (
            ["P12345\tIPR000001\tInterPro name\tPF00051\t10\t80"],
            """
<interpro id="IPR000001" type="Domain">
  <name>InterPro name</name>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00051" name=""/>
  </member_list>
</interpro>
""",
            "incomplete PFAM metadata",
        ),
        (
            [
                "P12345\tIPR000001\tInterPro name\tPF00051\t10\t80",
                "P12345\tIPR000003\tAnother InterPro name\tPF00051\t20\t90",
            ],
            """
<interpro id="IPR000001" type="Domain">
  <name>InterPro name</name>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00051" name="Kringle"/>
  </member_list>
</interpro>
<interpro id="IPR000003" type="Domain">
  <name>Another InterPro name</name>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00051" name="Conflicting name"/>
  </member_list>
</interpro>
""",
            "conflicting names",
        ),
        (
            ["P12345\tIPR000001\tInterPro name\tPF12\t10\t80"],
            None,
            "invalid PFAM member IDs",
        ),
        (
            ["P12345\tIPR999999\tMissing reference\tPF99999\t10\t80"],
            None,
            "incomplete PFAM metadata",
        ),
        (
            ["P12345\tIPR999999\tMismatched reference\tPF00051\t10\t80"],
            None,
            "incomplete PFAM metadata",
        ),
    ],
)
def test_build_pfam_tidy_rejects_invalid_contracts(
    tmp_path: Path,
    protein2ipr_rows: list[str],
    xml_entries: str | None,
    error: str,
) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(
        tmp_path,
        protein2ipr_rows=protein2ipr_rows,
        xml_entries=xml_entries,
    )

    with pytest.raises(ValueError, match=error):
        InterProDatabase.from_mapping_files(
            protein_to_interpro=file_protein2ipr,
            interpro_xml=file_xml,
        ).build_tidy(config="pfam")


def test_build_pfam_tidy_uses_xml_release_not_parent_directories(
    tmp_path: Path,
) -> None:
    file_protein2ipr, _file_xml_108 = write_interpro_snapshot(
        tmp_path / "snapshot-108",
        version="108.0",
    )
    _file_protein2ipr_109, file_xml = write_interpro_snapshot(
        tmp_path / "snapshot-109",
        version="109.0",
    )

    dataset = InterProDatabase.from_mapping_files(
        protein_to_interpro=file_protein2ipr,
        interpro_xml=file_xml,
    ).build_tidy(config="pfam")
    assert dataset.release_version == "109.0"
    assert dataset.release_version_source == "official_metadata"
    assert dataset.source_schema_version is None


def test_tidy_config_defaults_to_mapping(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=file_protein2ipr,
        interpro_xml=file_xml,
    )

    assert set(db.build_tidy().frames) == {"mapping"}
    assert set(db.build_tidy(config="mapping").frames) == {"mapping"}
    assert_frame_equal(
        db.build_tidy().frames["mapping"].collect(),
        db.build_tidy(config="mapping").frames["mapping"].collect(),
    )


def test_pfam_config_requires_xml(tmp_path: Path) -> None:
    file_protein2ipr, _file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDatabase.from_mapping_files(protein_to_interpro=file_protein2ipr)

    with pytest.raises(ValueError, match="XML file is required"):
        db.build_tidy(config="pfam")


def test_rejects_unknown_tidy_config(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDatabase.from_mapping_files(
        protein_to_interpro=file_protein2ipr,
        interpro_xml=file_xml,
    )

    with pytest.raises(ValueError, match="tidy config"):
        db.build_tidy(config="unknown")  # type: ignore[arg-type]


def test_standalone_pfam_api_is_not_exported() -> None:
    assert "build_pfam_tidy" not in interpro.__all__
    assert "write_pfam_tidy" not in interpro.__all__
