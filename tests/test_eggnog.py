from __future__ import annotations

import gzip
import shutil
import sqlite3
from pathlib import Path

import polars as pl
import pytest

from bioextract.eggnog import EggnogDb


def write_eggnog_fixture(tmp_path: Path) -> dict[str, Path]:
    file_db = tmp_path / "eggnog.db"
    with sqlite3.connect(file_db) as conn:
        conn.execute("create table prots (name text primary key, ogs text)")
        conn.execute(
            "create table og ("
            "og text, level text, description text, COG_categories text, "
            "primary key (og, level))"
        )
        conn.executemany(
            "insert into prots values (?, ?)",
            [
                ("9606.ENSP1", "OG0001@2759,OG0002@2759"),
                ("9606.ENSP2", "OG0003@2759"),
                ("9606.EMPTY", ""),
            ],
        )
        conn.executemany(
            "insert into og values (?, ?, ?, ?)",
            [
                ("OG0001", "2759", "alpha OG", "EG"),
                ("OG0002", "2759", "beta OG", "S"),
                ("OG0003", "2759", "gamma OG", "-"),
            ],
        )

    file_cog_fun = tmp_path / "cog-24.fun.tab"
    file_cog_fun.write_text(
        "\n".join(
            [
                "1\tINFORMATION STORAGE AND PROCESSING",
                "E\t3\tAAAAAA\tAmino acid transport and metabolism",
                "G\t3\tBBBBBB\tCarbohydrate transport and metabolism",
                "S\t4\tCCCCCC\tFunction unknown",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {"db": file_db, "cog_fun": file_cog_fun}


def test_extract_mapping_from_sqlite_expands_og_cog_categories(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggnogDb.from_files(
        file_eggnog_db=files["db"],
        file_cog_fun=files["cog_fun"],
    )

    df_mapping = db.extract_mapping()

    assert df_mapping.columns == [
        "EggnogProteinId",
        "EggnogOgId",
        "EggnogLevel",
        "CogCategory",
        "CogClass",
        "CogName",
        "OgDescription",
    ]
    assert df_mapping.to_dicts() == [
        {
            "EggnogProteinId": "9606.ENSP1",
            "EggnogOgId": "OG0001",
            "EggnogLevel": "2759",
            "CogCategory": "E",
            "CogClass": "3",
            "CogName": "Amino acid transport and metabolism",
            "OgDescription": "alpha OG",
        },
        {
            "EggnogProteinId": "9606.ENSP1",
            "EggnogOgId": "OG0001",
            "EggnogLevel": "2759",
            "CogCategory": "G",
            "CogClass": "3",
            "CogName": "Carbohydrate transport and metabolism",
            "OgDescription": "alpha OG",
        },
        {
            "EggnogProteinId": "9606.ENSP1",
            "EggnogOgId": "OG0002",
            "EggnogLevel": "2759",
            "CogCategory": "S",
            "CogClass": "4",
            "CogName": "Function unknown",
            "OgDescription": "beta OG",
        },
    ]


def test_select_ids_queries_subset_and_reports_unmapped(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggnogDb.from_files(
        file_eggnog_db=files["db"],
        file_cog_fun=files["cog_fun"],
    )

    selection = db.select_ids(
        ["9606.ENSP1", "9606.MISSING"],
        kind_input_id="eggnog_protein",
    )

    df_mapping = selection.extract_mapping()
    assert df_mapping.select("InputId", "KindInputId", "CogCategory").to_dicts() == [
        {
            "InputId": "9606.ENSP1",
            "KindInputId": "eggnog_protein",
            "CogCategory": "E",
        },
        {
            "InputId": "9606.ENSP1",
            "KindInputId": "eggnog_protein",
            "CogCategory": "G",
        },
        {
            "InputId": "9606.ENSP1",
            "KindInputId": "eggnog_protein",
            "CogCategory": "S",
        },
    ]
    assert selection.extract_unmapped_input_ids().to_dicts() == [
        {"InputId": "9606.MISSING"}
    ]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggnogDb.from_files(
        file_eggnog_db=files["db"],
        file_cog_fun=files["cog_fun"],
    )

    selection = db.select_groups(
        {"up": ["9606.ENSP1"], "down": ["9606.MISSING"]},
        kind_input_id="eggnog_protein",
    )

    assert selection.extract_mapping().columns[:3] == [
        "GroupId",
        "InputId",
        "KindInputId",
    ]
    assert selection.extract_mapping().select("GroupId", "CogCategory").to_dicts() == [
        {"GroupId": "up", "CogCategory": "E"},
        {"GroupId": "up", "CogCategory": "G"},
        {"GroupId": "up", "CogCategory": "S"},
    ]
    assert selection.extract_unmapped_input_ids().to_dicts() == [
        {"GroupId": "down", "InputId": "9606.MISSING"}
    ]


def test_build_tidy_writes_mapping_parquet_and_manifest(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggnogDb.from_files(
        file_eggnog_db=files["db"],
        file_cog_fun=files["cog_fun"],
    )

    report = db.write_tidy(tmp_path / "out", should_write_manifest=True)

    assert report.manifest is not None
    assert report.manifest["schema_version"] == "eggnog-mapping-v0.1"
    assert [asset.path for asset in report.assets] == ["mapping.parquet"]
    assert pl.read_parquet(tmp_path / "out" / "mapping.parquet").height == 3


def test_gzip_sqlite_is_decompressed_to_tmp_dir(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    file_gz = tmp_path / "eggnog.db.gz"
    with files["db"].open("rb") as handle_in:
        with gzip.open(file_gz, "wb") as handle_out:
            shutil.copyfileobj(handle_in, handle_out)

    db = EggnogDb.from_files(
        file_eggnog_db=file_gz,
        file_cog_fun=files["cog_fun"],
        dir_tmp=tmp_path / "tmp",
    )

    assert (
        db.select_ids(
            ["9606.ENSP1"],
            kind_input_id="eggnog_protein",
        )
        .extract_mapping()
        .height
        == 3
    )


def test_validates_kind_input_id(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggnogDb.from_files(file_eggnog_db=files["db"])

    with pytest.raises(ValueError, match="kind_input_id"):
        db.select_ids(["9606.ENSP1"], kind_input_id="uniprot")  # type: ignore[arg-type]
