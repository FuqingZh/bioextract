from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from bioextract.uniprot import UniProtDatabase
from bioextract.uniprot.constant import (
    COLS_IDMAPPING_SELECTED,
    COLS_SUBCELLULAR_LOCATION,
    SCHEMA_SUBCELLULAR_LOCATION,
    SCHEMA_VERSION_SUBCELLULAR_LOCATION,
)
from bioextract.uniprot.util import read_subcellular_location_frame


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


def write_uniprot_subcellular_dat_fixture(
    tmp_path: Path,
    *,
    should_gzip: bool = True,
) -> Path:
    suffix = ".dat.gz" if should_gzip else ".dat"
    file_in = tmp_path / f"uniprot_sprot_subcellular{suffix}"
    text = """ID   TEST1_HUMAN              Reviewed;         100 AA.
AC   P12345; Q11111;
DE   RecName: Full=Cellular tumor antigen p53;
GN   Name=TP53; Synonyms=P53;
CC   -!- SUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:123456}.
CC       Nucleus {ECO:0000305}. Note=Shuttles between cytoplasm and nucleus.
CC   -!- FUNCTION: DNA-binding transcription factor.
DR   eggNOG; KOG0001; Eukaryota.
//
ID   TEST2_HUMAN              Reviewed;         200 AA.
AC   Q9Y243;
DE   RecName: Full=RAC-gamma serine/threonine-protein kinase {ECO:0000303|Ref.6};
GN   ORFNames=CG14996;
CC   -!- SUBCELLULAR LOCATION: Membrane.
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
        UniProtDatabase.from_files(id_mapping=file_in)
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

    df_xref = UniProtDatabase.from_dat(
        path=file_dat,
        source_database="sprot",
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
    db = UniProtDatabase.from_dat(path=file_dat, source_database="sprot")

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


def test_subcellular_location_schema_constants_are_declared() -> None:
    assert SCHEMA_VERSION_SUBCELLULAR_LOCATION == "uniprot-subcellular-location-v0.1"
    assert COLS_SUBCELLULAR_LOCATION == [
        "UniProtId",
        "PrimaryUniProtId",
        "UniProtEntryName",
        "GeneName",
        "ProteinName",
        "SubcellularLocation",
        "SubcellularLocationNote",
        "EvidenceCode",
        "EvidenceSource",
        "EvidenceId",
        "SourceDb",
    ]
    assert list(SCHEMA_SUBCELLULAR_LOCATION) == COLS_SUBCELLULAR_LOCATION


def test_read_subcellular_location_frame_from_dat_gzip(tmp_path: Path) -> None:
    file_dat = write_uniprot_subcellular_dat_fixture(tmp_path)

    df_subcellular = read_subcellular_location_frame(file_dat, source_db="sprot")

    assert df_subcellular.columns == COLS_SUBCELLULAR_LOCATION
    assert df_subcellular.to_dicts() == [
        {
            "UniProtId": "P12345",
            "PrimaryUniProtId": "P12345",
            "UniProtEntryName": "TEST1_HUMAN",
            "GeneName": "TP53",
            "ProteinName": "Cellular tumor antigen p53",
            "SubcellularLocation": "Cytoplasm",
            "SubcellularLocationNote": "Shuttles between cytoplasm and nucleus",
            "EvidenceCode": "ECO:0000269",
            "EvidenceSource": "PubMed",
            "EvidenceId": "123456",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "P12345",
            "PrimaryUniProtId": "P12345",
            "UniProtEntryName": "TEST1_HUMAN",
            "GeneName": "TP53",
            "ProteinName": "Cellular tumor antigen p53",
            "SubcellularLocation": "Nucleus",
            "SubcellularLocationNote": "Shuttles between cytoplasm and nucleus",
            "EvidenceCode": "ECO:0000305",
            "EvidenceSource": None,
            "EvidenceId": None,
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "UniProtEntryName": "TEST1_HUMAN",
            "GeneName": "TP53",
            "ProteinName": "Cellular tumor antigen p53",
            "SubcellularLocation": "Cytoplasm",
            "SubcellularLocationNote": "Shuttles between cytoplasm and nucleus",
            "EvidenceCode": "ECO:0000269",
            "EvidenceSource": "PubMed",
            "EvidenceId": "123456",
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q11111",
            "PrimaryUniProtId": "P12345",
            "UniProtEntryName": "TEST1_HUMAN",
            "GeneName": "TP53",
            "ProteinName": "Cellular tumor antigen p53",
            "SubcellularLocation": "Nucleus",
            "SubcellularLocationNote": "Shuttles between cytoplasm and nucleus",
            "EvidenceCode": "ECO:0000305",
            "EvidenceSource": None,
            "EvidenceId": None,
            "SourceDb": "sprot",
        },
        {
            "UniProtId": "Q9Y243",
            "PrimaryUniProtId": "Q9Y243",
            "UniProtEntryName": "TEST2_HUMAN",
            "GeneName": None,
            "ProteinName": "RAC-gamma serine/threonine-protein kinase",
            "SubcellularLocation": "Membrane",
            "SubcellularLocationNote": None,
            "EvidenceCode": None,
            "EvidenceSource": None,
            "EvidenceId": None,
            "SourceDb": "sprot",
        },
    ]


def test_read_subcellular_location_frame_from_plain_dat(tmp_path: Path) -> None:
    file_dat = write_uniprot_subcellular_dat_fixture(tmp_path, should_gzip=False)

    df_subcellular = read_subcellular_location_frame(file_dat, source_db="sprot")

    assert df_subcellular.height == 5
    assert df_subcellular.select("UniProtId", "SubcellularLocation").to_dicts()[-1] == {
        "UniProtId": "Q9Y243",
        "SubcellularLocation": "Membrane",
    }


def test_extract_subcellular_location_from_dat_handle(tmp_path: Path) -> None:
    file_dat = write_uniprot_subcellular_dat_fixture(tmp_path)

    df_subcellular = UniProtDatabase.from_dat(
        path=file_dat,
        source_database="sprot",
    ).extract_subcellular_location()

    assert df_subcellular.columns == COLS_SUBCELLULAR_LOCATION
    assert df_subcellular.height == 5
    assert sorted(
        df_subcellular.select("PrimaryUniProtId").unique().to_series().to_list()
    ) == ["P12345", "Q9Y243"]


def test_write_subcellular_location_parquet_writes_data_without_sidecar(
    tmp_path: Path,
) -> None:
    file_dat = write_uniprot_subcellular_dat_fixture(tmp_path)
    db = UniProtDatabase.from_dat(path=file_dat, source_database="sprot")

    path = tmp_path / "subcellular_location.parquet"
    result = db.write_subcellular_location_parquet(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    df_written = pl.read_parquet(path)
    assert df_written.columns == [
        "uniprot_id",
        "primary_uniprot_id",
        "uniprot_entry_name",
        "gene_name",
        "protein_name",
        "subcellular_location",
        "subcellular_location_note",
        "evidence_code",
        "evidence_source",
        "evidence_id",
        "source_db",
    ]
    assert df_written.height == 5


def test_write_eggnog_xref_parquet_writes_mapping_without_sidecar(
    tmp_path: Path,
) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path)
    db = UniProtDatabase.from_dat(path=file_dat, source_database="sprot")

    path = tmp_path / "eggnog_xref.parquet"
    result = db.write_eggnog_xref_parquet(path)

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()
    assert pl.read_parquet(path).height == 5


def test_uniprot_dat_and_idmapping_methods_are_separated(tmp_path: Path) -> None:
    file_dat = write_uniprot_dat_fixture(tmp_path)
    db_dat = UniProtDatabase.from_dat(path=file_dat, source_database="sprot")

    with pytest.raises(ValueError, match="idmapping"):
        db_dat.extract_mapping()

    file_idmapping = write_idmapping_fixture(tmp_path)
    db_idmapping = UniProtDatabase.from_files(id_mapping=file_idmapping)
    with pytest.raises(ValueError, match="eggNOG"):
        db_idmapping.extract_eggnog_xref()
    with pytest.raises(ValueError, match="subcellular"):
        db_idmapping.extract_subcellular_location()


def test_write_parquet_writes_single_file_by_default(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    path = tmp_path / "uniprot.parquet"

    result = (
        UniProtDatabase.from_files(id_mapping=file_in)
        .with_taxids("9606", "10090")
        .write_parquet(path)
    )

    assert result.path == path
    assert not (tmp_path / "manifest.json").exists()

    df_hsa = (
        UniProtDatabase.from_files(id_mapping=path)
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
        UniProtDatabase.from_files(id_mapping=file_in)
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
        UniProtDatabase.from_files(id_mapping=dir_out)
        .with_taxids("9606")
        .extract_mapping()
    )

    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_parquet_all_requires_explicit_flag(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    db = UniProtDatabase.from_files(id_mapping=file_in)

    with pytest.raises(ValueError, match="allow_all_taxa=True"):
        db.write_parquet(tmp_path / "blocked.parquet")

    result = db.write_parquet(tmp_path / "all.parquet", allow_all_taxa=True)
    assert result.path == tmp_path / "all.parquet"


def test_single_taxid_can_write_single_parquet_and_read_it(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path, should_gzip=False)
    path = tmp_path / "single.parquet"

    UniProtDatabase.from_files(id_mapping=file_in).with_taxids("9606").write_parquet(
        path
    )

    df_mapping = UniProtDatabase.from_files(id_mapping=path).extract_mapping()
    assert df_mapping.height == 2
    assert df_mapping["TaxId"].unique().to_list() == ["9606"]


def test_write_parquet_preserves_existing_output_by_default(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)
    path = tmp_path / "uniprot.parquet"
    path.write_bytes(b"existing")
    db = UniProtDatabase.from_files(id_mapping=file_in).with_taxids("9606")

    with pytest.raises(FileExistsError):
        db.write_parquet(path)
    assert path.read_bytes() == b"existing"


def test_validate_schema_rejects_bad_parquet(tmp_path: Path) -> None:
    file_bad = tmp_path / "bad.parquet"
    pl.DataFrame({"UniProtId": ["P04637"]}).write_parquet(file_bad)

    with pytest.raises(ValueError, match="missing required columns"):
        UniProtDatabase.from_files(id_mapping=file_bad).validate_schema()


def test_from_files_rejects_missing_unsupported_and_empty_hive(
    tmp_path: Path,
) -> None:
    with pytest.raises(FileNotFoundError):
        UniProtDatabase.from_files(id_mapping=tmp_path / "missing.tab.gz")

    file_unsupported = tmp_path / "mapping.txt"
    file_unsupported.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported UniProt"):
        UniProtDatabase.from_files(id_mapping=file_unsupported)

    dir_empty = tmp_path / "empty"
    dir_empty.mkdir()
    with pytest.raises(ValueError, match="contains no parquet"):
        UniProtDatabase.from_files(id_mapping=dir_empty)

    dir_no_parquet = tmp_path / "no_parquet"
    (dir_no_parquet / "TaxId=9606").mkdir(parents=True)
    (dir_no_parquet / "manifest.json").write_text("{}", encoding="utf-8")
    (dir_no_parquet / "TaxId=9606" / "note.txt").write_text("x", encoding="utf-8")
    with pytest.raises(ValueError, match="contains no parquet"):
        UniProtDatabase.from_files(id_mapping=dir_no_parquet)


def test_taxid_validation(tmp_path: Path) -> None:
    file_in = write_idmapping_fixture(tmp_path)

    with pytest.raises(ValueError, match="TaxId values"):
        UniProtDatabase.from_files(id_mapping=file_in).with_taxids("9606", "")
