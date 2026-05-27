from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from bioextract.uniprot import UniprotDb, UniprotResourceLimits
from bioextract.uniprot.constant import COLS_IDMAPPING_SELECTED


def write_idmapping_fixture(tmp_path: Path, *, should_gzip: bool = True) -> Path:
    suffix = ".tab.gz" if should_gzip else ".tab"
    file_in = tmp_path / f"idmapping_selected{suffix}"
    rows = [
        [
            "P04637",
            "P53_HUMAN",
            "7157",
            "NP_000537.3",
            "",
            "1TUP",
            "GO:0006915; GO:0003677",
            "UniRef100_P04637",
            "UniRef90_P04637",
            "UniRef50_P04637",
            "UPI000002ED67",
            "",
            "9606",
            "",
            "",
            "123456",
            "X54156",
            "CAA38295.1",
            "ENSG00000141510",
            "ENST00000269305",
            "ENSP00000269305",
            "",
        ],
        [
            "Q9Y243",
            "AKT3_HUMAN",
            "10000",
            "NP_005456.1",
            "",
            "",
            "GO:0004672",
            "UniRef100_Q9Y243",
            "UniRef90_Q9Y243",
            "UniRef50_Q9Y243",
            "UPI000013D9F8",
            "",
            "9606",
            "",
            "",
            "",
            "",
            "",
            "ENSG00000117020",
            "",
            "",
            "",
        ],
        [
            "P31750",
            "AKT1_MOUSE",
            "11651",
            "NP_033782.1",
            "",
            "",
            "GO:0004672",
            "UniRef100_P31750",
            "UniRef90_P31750",
            "UniRef50_P31750",
            "UPI0000000001",
            "",
            "10090",
            "",
            "",
            "",
            "",
            "",
            "ENSMUSG00000001729",
            "",
            "",
            "",
        ],
    ]
    text = "\n".join("\t".join(row) for row in rows) + "\n"
    if should_gzip:
        with gzip.open(file_in, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        file_in.write_text(text, encoding="utf-8")
    return file_in


def test_extract_mapping_filters_taxids_from_raw_gzip(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    df_mapping = (
        UniprotDb.from_files(file_idmapping_selected=file_in)
        .with_taxids("9606")
        .extract_mapping()
    )

    assert df_mapping.columns == COLS_IDMAPPING_SELECTED
    assert df_mapping.select("UniProtId", "GeneId", "TaxId").to_dicts() == [
        {"UniProtId": "P04637", "GeneId": "7157", "TaxId": "9606"},
        {"UniProtId": "Q9Y243", "GeneId": "10000", "TaxId": "9606"},
    ]


def test_write_tidy_defaults_to_hive_and_can_read_hive_dataset(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"

    report = (
        UniprotDb.from_files(file_idmapping_selected=file_in)
        .with_taxids("9606", "10090")
        .write_tidy(dir_out, should_write_manifest=True)
    )

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "uniprot-idmapping-selected-v0.1"
    assert sorted(
        path.relative_to(dir_out).as_posix() for path in dir_out.rglob("*.parquet")
    ) == [
        "TaxId=10090/00000000.parquet",
        "TaxId=9606/00000000.parquet",
    ]

    df_hsa = (
        UniprotDb.from_files(file_idmapping_selected=dir_out)
        .with_taxids("9606")
        .extract_mapping()
    )
    assert df_hsa["TaxId"].unique().to_list() == ["9606"]
    assert df_hsa.height == 2


def test_hive_dataset_allows_manifest_and_nested_non_parquet_files(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"
    (
        UniprotDb.from_files(file_idmapping_selected=file_in)
        .with_taxids("9606")
        .write_tidy(dir_out, should_write_manifest=True)
    )
    (dir_out / "README.txt").write_text("metadata", encoding="utf-8")
    (dir_out / "TaxId=9606" / "note.txt").write_text("metadata", encoding="utf-8")

    df_mapping = (
        UniprotDb.from_files(file_idmapping_selected=dir_out)
        .with_taxids("9606")
        .extract_mapping()
    )

    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_tidy_all_requires_explicit_flag_and_writes_hive(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    db = UniprotDb.from_files(file_idmapping_selected=file_in)

    with pytest.raises(ValueError, match="should_allow_all=True"):
        db.write_tidy(tmp_path / "blocked")

    report = db.write_tidy(tmp_path / "all", should_allow_all=True)

    assert sorted(asset["path"] for asset in report.assets) == [
        "TaxId=10090/00000000.parquet",
        "TaxId=9606/00000000.parquet",
    ]


def test_single_taxid_can_write_single_parquet_and_read_it(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path, should_gzip=False)
    file_out = tmp_path / "single"

    UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
        "9606"
    ).write_tidy(file_out, should_write_hive=False)

    df_mapping = UniprotDb.from_files(
        file_idmapping_selected=file_out / "mapping.parquet"
    ).extract_mapping()
    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_tidy_rejects_non_empty_output_unless_overwrite(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"
    dir_out.mkdir()
    (dir_out / "old.txt").write_text("old", encoding="utf-8")
    db = UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids("9606")

    with pytest.raises(FileExistsError, match="not empty"):
        db.write_tidy(dir_out)

    db.write_tidy(dir_out, should_overwrite=True)

    assert not (dir_out / "old.txt").exists()
    assert sorted(
        path.relative_to(dir_out).as_posix() for path in dir_out.rglob("*.parquet")
    ) == [
        "TaxId=9606/00000000.parquet",
    ]


def test_validate_schema_rejects_bad_parquet(tmp_path: Path) -> None:
    file_bad = tmp_path / "bad.parquet"
    pl.DataFrame({"UniProtId": ["P04637"]}).write_parquet(file_bad)

    with pytest.raises(ValueError, match="missing required columns"):
        UniprotDb.from_files(file_idmapping_selected=file_bad).validate_schema()


def test_from_files_rejects_missing_unsupported_and_empty_hive(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        UniprotDb.from_files(file_idmapping_selected=tmp_path / "missing.tab.gz")

    file_unsupported = tmp_path / "mapping.txt"
    file_unsupported.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported UniProt"):
        UniprotDb.from_files(file_idmapping_selected=file_unsupported)

    dir_empty = tmp_path / "empty"
    dir_empty.mkdir()
    with pytest.raises(ValueError, match="contains no parquet"):
        UniprotDb.from_files(file_idmapping_selected=dir_empty)

    dir_no_parquet = tmp_path / "no_parquet"
    (dir_no_parquet / "TaxId=9606").mkdir(parents=True)
    (dir_no_parquet / "manifest.json").write_text("{}", encoding="utf-8")
    (dir_no_parquet / "TaxId=9606" / "note.txt").write_text("x", encoding="utf-8")
    with pytest.raises(ValueError, match="contains no parquet"):
        UniprotDb.from_files(file_idmapping_selected=dir_no_parquet)


def test_resource_limits_and_taxid_validation(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    with pytest.raises(ValueError, match="size limit"):
        UniprotDb.from_files(
            file_idmapping_selected=file_in,
            limits=UniprotResourceLimits(file_idmapping_selected_bytes_max=1),
        )

    with pytest.raises(ValueError, match="TaxId values"):
        UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids("9606", "")
