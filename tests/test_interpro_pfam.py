from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

import bioextract.interpro as interpro
from bioextract.interpro import InterProDb


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


def test_write_pfam_tidy_emits_compact_assets_and_hashes(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    dir_out = tmp_path / "out"

    db = InterProDb.from_mapping_files(
        file_protein2ipr=file_protein2ipr,
        file_interpro_xml=file_xml,
    )
    report = db.write_tidy(
        dir_out,
        config="pfam",
        should_write_manifest=True,
        should_hash_sources=True,
        should_hash_assets=True,
    )

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "interpro-pfam-v0.1"
    assert [asset.path for asset in report.assets] == [
        "protein_term.parquet",
        "term.parquet",
        "term_xref.parquet",
    ]
    assert all(
        isinstance(sha256 := source.get("sha256"), str) and len(sha256) == 64
        for source in report.manifest["sources"]
    )
    assert all(
        isinstance(asset["sha256"], str) and len(asset["sha256"]) == 64
        for asset in report.manifest["assets"]
    )

    df_protein_term = pl.read_parquet(dir_out / "protein_term.parquet").sort(
        "UniProtId", "PfamId"
    )
    assert df_protein_term.schema == {
        "UniProtId": pl.String,
        "PfamId": pl.String,
    }
    assert df_protein_term.to_dicts() == [
        {"UniProtId": "P12345", "PfamId": "PF00051"},
        {"UniProtId": "P12345", "PfamId": "PF00069"},
        {"UniProtId": "Q9Y243", "PfamId": "PF00069"},
    ]

    df_term = pl.read_parquet(dir_out / "term.parquet")
    assert df_term.to_dicts() == [
        {"PfamId": "PF00051", "PfamName": "Kringle"},
        {"PfamId": "PF00069", "PfamName": "Pkinase"},
    ]
    assert "InterPro Kringle domain" not in df_term.get_column("PfamName").to_list()

    df_term_xref = pl.read_parquet(dir_out / "term_xref.parquet")
    assert df_term_xref.to_dicts() == [
        {
            "PfamId": "PF00051",
            "InterProId": "IPR000001",
            "InterProName": "InterPro Kringle domain",
            "InterProType": "Domain",
        },
        {
            "PfamId": "PF00069",
            "InterProId": "IPR000002",
            "InterProName": "Protein kinase domain",
            "InterProType": "Homologous_superfamily",
        },
    ]


def test_build_pfam_tidy_omits_source_hashes_by_default(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)

    db = InterProDb.from_mapping_files(
        file_protein2ipr=file_protein2ipr,
        file_interpro_xml=file_xml,
    )
    dataset = db.build_tidy(config="pfam")
    manifest = dataset.build_manifest([])

    assert dataset.schema_version == "interpro-pfam-v0.1"
    assert all("sha256" not in source for source in manifest["sources"])
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
        InterProDb.from_mapping_files(
            file_protein2ipr=file_protein2ipr,
            file_interpro_xml=file_xml,
        ).build_tidy(config="pfam")


def test_build_pfam_tidy_rejects_cross_version_inputs(tmp_path: Path) -> None:
    file_protein2ipr, _file_xml_108 = write_interpro_snapshot(
        tmp_path / "snapshot-108",
        version="108.0",
    )
    _file_protein2ipr_109, file_xml = write_interpro_snapshot(
        tmp_path / "snapshot-109",
        version="109.0",
    )

    with pytest.raises(ValueError, match="same snapshot directory"):
        InterProDb.from_mapping_files(
            file_protein2ipr=file_protein2ipr,
            file_interpro_xml=file_xml,
        ).build_tidy(config="pfam")


def test_tidy_config_defaults_to_mapping(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=file_protein2ipr,
        file_interpro_xml=file_xml,
    )

    assert set(db.build_tidy().frames) == {"mapping"}
    assert set(db.build_tidy(config="mapping").frames) == {"mapping"}
    assert_frame_equal(
        db.build_tidy().frames["mapping"].collect(),
        db.build_tidy(config="mapping").frames["mapping"].collect(),
    )


def test_pfam_config_requires_xml(tmp_path: Path) -> None:
    file_protein2ipr, _file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDb.from_mapping_files(file_protein2ipr=file_protein2ipr)

    with pytest.raises(ValueError, match="XML file is required"):
        db.build_tidy(config="pfam")


def test_rejects_unknown_tidy_config(tmp_path: Path) -> None:
    file_protein2ipr, file_xml = write_interpro_snapshot(tmp_path)
    db = InterProDb.from_mapping_files(
        file_protein2ipr=file_protein2ipr,
        file_interpro_xml=file_xml,
    )

    with pytest.raises(ValueError, match="tidy config"):
        db.build_tidy(config="unknown")  # type: ignore[arg-type]


def test_standalone_pfam_api_is_not_exported() -> None:
    assert "build_pfam_tidy" not in interpro.__all__
    assert "write_pfam_tidy" not in interpro.__all__
