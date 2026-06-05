from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

import bioextract.uniprot.uniprot as uniprot_module
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


def write_uniprot_dat_fixture(tmp_path: Path, *, should_gzip: bool = True) -> Path:
    suffix = ".dat.gz" if should_gzip else ".dat"
    file_in = tmp_path / f"uniprot_sprot{suffix}"
    text = """ID   TEST1_HUMAN              Reviewed;         100 AA.
AC   P12345; Q11111;
DR   eggNOG; KOG0001; Eukaryota.
DR   eggNOG; ENOG502ABC; Metazoa.
//
ID   TEST2_BACTERIA           Reviewed;         200 AA.
AC   Q9Y243;
DR   eggNOG; COG1028; Bacteria.
//
ID   TEST3_HUMAN              Reviewed;         300 AA.
AC   P31750;
//
"""
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


def test_extract_eggnog_xref_from_dat_gzip(tmp_path: Path) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path)

    df_xref = UniprotDb.from_dat(
        file_dat=file_dat,
        source_db="sprot",
    ).extract_eggnog_xref()

    assert df_xref.to_dicts() == [
        {
            "UniProtId": "P12345",
            "PrimaryUniProtId": "P12345",
            "IsPrimaryAccession": True,
            "EggnogOgId": "ENOG502ABC",
            "EggnogLevel": "Metazoa",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "P12345",
            "PrimaryUniProtId": "P12345",
            "IsPrimaryAccession": True,
            "EggnogOgId": "KOG0001",
            "EggnogLevel": "Eukaryota",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "IsPrimaryAccession": False,
            "EggnogOgId": "ENOG502ABC",
            "EggnogLevel": "Metazoa",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "IsPrimaryAccession": False,
            "EggnogOgId": "KOG0001",
            "EggnogLevel": "Eukaryota",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q9Y243",
            "PrimaryUniProtId": "Q9Y243",
            "IsPrimaryAccession": True,
            "EggnogOgId": "COG1028",
            "EggnogLevel": "Bacteria",
            "SourceDb": "sprot",
        },
    ]


def test_select_eggnog_xref_ids_filters_accessions(tmp_path: Path) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path, should_gzip=False)
    db = UniprotDb.from_dat(file_dat=file_dat, source_db="sprot")

    df_xref = db.select_eggnog_xref_ids(["Q11111", "MISSING"])

    assert df_xref.select("UniProtId", "PrimaryUniProtId", "EggnogOgId").to_dicts() == [
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "EggnogOgId": "ENOG502ABC",
        },
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "EggnogOgId": "KOG0001",
        },
    ]


def test_write_eggnog_xref_tidy_writes_mapping_parquet_and_manifest(
    tmp_path: Path,
) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path)
    db = UniprotDb.from_dat(file_dat=file_dat, source_db="sprot")

    report = db.write_eggnog_xref_tidy(
        tmp_path / "eggnog-xref",
        should_write_manifest=True,
    )

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "uniprot-eggnog-xref-v0.1"
    assert [asset["path"] for asset in report.assets] == ["mapping.parquet"]
    assert pl.read_parquet(tmp_path / "eggnog-xref" / "mapping.parquet").height == 5


def test_uniprot_dat_and_idmapping_methods_are_separated(tmp_path: Path) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path)
    db_dat = UniprotDb.from_dat(file_dat=file_dat, source_db="sprot")

    with pytest.raises(ValueError, match="idmapping"):
        db_dat.extract_mapping()

    file_idmapping = write_idmapping_fixture(tmp_path)
    db_idmapping = UniprotDb.from_files(file_idmapping_selected=file_idmapping)
    with pytest.raises(ValueError, match="eggNOG"):
        db_idmapping.extract_eggnog_xref()


def test_write_tidy_writes_single_parquet_by_default(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"

    report = (
        UniprotDb.from_files(file_idmapping_selected=file_in)
        .with_taxids("9606", "10090")
        .write_tidy(dir_out, should_write_manifest=True)
    )

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "uniprot-idmapping-selected-v0.1"
    assert report.assets == (
        {
            "path": "mapping.parquet",
            "kind": "canonical",
            "row_count": None,
            "is_optional": False,
        },
    )
    assert (dir_out / "mapping.parquet").is_file()

    df_hsa = (
        UniprotDb.from_files(file_idmapping_selected=dir_out / "mapping.parquet")
        .with_taxids("9606")
        .extract_mapping()
    )
    assert df_hsa["TaxId"].unique().to_list() == ["9606"]
    assert df_hsa.height == 2


def test_hive_dataset_allows_manifest_and_nested_non_parquet_files(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    df_fixture = (
        UniprotDb.from_files(file_idmapping_selected=file_in)
        .with_taxids("9606")
        .extract_mapping()
    )
    dir_out = tmp_path / "tidy"
    dir_taxid = dir_out / "TaxId=9606"
    dir_taxid.mkdir(parents=True)
    df_fixture.write_parquet(dir_taxid / "00000000.parquet")
    (dir_out / "README.txt").write_text("metadata", encoding="utf-8")
    (dir_taxid / "note.txt").write_text("metadata", encoding="utf-8")

    df_mapping = (
        UniprotDb.from_files(file_idmapping_selected=dir_out)
        .with_taxids("9606")
        .extract_mapping()
    )

    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_tidy_all_requires_explicit_flag_and_writes_single_parquet(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    db = UniprotDb.from_files(file_idmapping_selected=file_in)

    with pytest.raises(ValueError, match="should_allow_all=True"):
        db.write_tidy(tmp_path / "blocked")

    report = db.write_tidy(tmp_path / "all", should_allow_all=True)

    assert report.assets == (
        {
            "path": "mapping.parquet",
            "kind": "canonical",
            "row_count": None,
            "is_optional": False,
        },
    )
    assert (tmp_path / "all" / "mapping.parquet").is_file()


def test_single_taxid_can_write_single_parquet_and_read_it(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path, should_gzip=False)
    file_out = tmp_path / "single"

    UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
        "9606"
    ).write_tidy(file_out)

    df_mapping = UniprotDb.from_files(
        file_idmapping_selected=file_out / "mapping.parquet"
    ).extract_mapping()
    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_tidy_applies_existing_output_policy(
    tmp_path: Path,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"
    dir_out.mkdir()
    (dir_out / "old.txt").write_text("old", encoding="utf-8")
    db = UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids("9606")

    with pytest.raises(FileExistsError, match="not empty"):
        db.write_tidy(dir_out)

    report_skip = db.write_tidy(dir_out, policy_existing="skip")

    assert (dir_out / "old.txt").read_text(encoding="utf-8") == "old"
    assert report_skip.assets == (
        {
            "path": "mapping.parquet",
            "kind": "canonical",
            "row_count": None,
            "is_optional": False,
        },
    )

    db.write_tidy(dir_out, policy_existing="overwrite")

    assert not (dir_out / "old.txt").exists()
    assert sorted(
        path.relative_to(dir_out).as_posix() for path in dir_out.rglob("*.parquet")
    ) == [
        "mapping.parquet",
    ]


def test_write_tidy_rejects_invalid_existing_output_policy(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    with pytest.raises(ValueError, match="policy_existing"):
        UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
            "9606"
        ).write_tidy(tmp_path / "tidy", policy_existing="replace")  # type: ignore[arg-type]


def test_write_tidy_skip_can_return_existing_manifest(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_out = tmp_path / "tidy"
    db = UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids("9606")

    db.write_tidy(dir_out, should_write_manifest=True)
    report_skip = db.write_tidy(
        dir_out,
        policy_existing="skip",
        should_write_manifest=True,
    )

    assert report_skip.manifest is not None
    assert report_skip.manifest["schema_version"] == "uniprot-idmapping-selected-v0.1"
    assert report_skip.assets == tuple(report_skip.manifest["assets"])


def test_write_tidy_accepts_zstd_compression_level(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
        "9606"
    ).write_tidy(tmp_path / "tidy", level_compression=1)

    assert (tmp_path / "tidy" / "mapping.parquet").is_file()


def test_write_tidy_rejects_all_taxid_ceph_without_local_tmp(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    db = UniprotDb.from_files(file_idmapping_selected=file_in)

    with pytest.raises(ValueError, match="explicit local dir_tmp"):
        db.write_tidy("/cephfs_data/example/uniprot/tidy", should_allow_all=True)

    with pytest.raises(ValueError, match="must not be under /cephfs_data"):
        db.write_tidy(
            "/cephfs_data/example/uniprot/tidy",
            should_allow_all=True,
            dir_tmp="/cephfs_data/tmp",
        )


def test_write_tidy_uses_explicit_tmp_and_publishes_output(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    dir_tmp = tmp_path / "scratch"
    dir_out = tmp_path / "tidy"

    UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
        "9606"
    ).write_tidy(dir_out, dir_tmp=dir_tmp)

    assert (dir_out / "mapping.parquet").is_file()
    assert sorted(path.name for path in dir_tmp.iterdir()) == []


def test_write_tidy_resource_monitor_stops_on_rss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    monkeypatch.setattr(
        uniprot_module,
        "_sample_process_resources",
        lambda: uniprot_module._ResourceSample(
            size_rss_mb=2,
            size_threads=1,
            count_d_state_threads=0,
        ),
    )

    with pytest.raises(RuntimeError, match="RSS stop threshold"):
        UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
            "9606"
        ).write_tidy(tmp_path / "tidy", size_rss_stop_gb=0.001)


def test_write_tidy_resource_monitor_stops_on_threads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    monkeypatch.setattr(
        uniprot_module,
        "_sample_process_resources",
        lambda: uniprot_module._ResourceSample(
            size_rss_mb=1,
            size_threads=3,
            count_d_state_threads=0,
        ),
    )

    with pytest.raises(RuntimeError, match="thread stop threshold"):
        UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
            "9606"
        ).write_tidy(tmp_path / "tidy", size_threads_stop=2)


def test_write_tidy_resource_monitor_stops_on_repeated_d_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    monkeypatch.setattr(
        uniprot_module,
        "_sample_process_resources",
        lambda: uniprot_module._ResourceSample(
            size_rss_mb=1,
            size_threads=1,
            count_d_state_threads=1,
        ),
    )

    with pytest.raises(RuntimeError, match="D-state"):
        UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
            "9606"
        ).write_tidy(tmp_path / "tidy", count_d_state_stop=2)


def test_write_tidy_can_disable_resource_monitor(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    monkeypatch.setattr(
        uniprot_module,
        "_sample_process_resources",
        lambda: uniprot_module._ResourceSample(
            size_rss_mb=10_000,
            size_threads=10_000,
            count_d_state_threads=10_000,
        ),
    )

    UniprotDb.from_files(file_idmapping_selected=file_in).with_taxids(
        "9606"
    ).write_tidy(tmp_path / "tidy", should_monitor_resources=False)

    assert (tmp_path / "tidy" / "mapping.parquet").is_file()


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
