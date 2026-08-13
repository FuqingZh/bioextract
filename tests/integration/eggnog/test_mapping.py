from __future__ import annotations

import gzip
import shutil
import sqlite3
import warnings
from collections.abc import Iterable, Iterator
from pathlib import Path

import polars as pl
import pytest

from bioextract.eggnog import EggNOGDatabase
from bioextract.eggnog import util as eggnog_util


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


def test_selected_mapping_from_sqlite_expands_og_cog_categories(
    tmp_path: Path,
) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggNOGDatabase.from_sqlite(
        files["db"],
        cog_functions=files["cog_fun"],
    )

    df_mapping = (
        db.select_ids(["9606.ENSP1"])
        .mappings()
        .collect()
        .select(
            "name",
            "og",
            "level",
            "description",
            "COG_categories",
            "cog_category",
            "cog_class",
            "cog_name",
        )
    )

    assert df_mapping.columns == [
        "name",
        "og",
        "level",
        "description",
        "COG_categories",
        "cog_category",
        "cog_class",
        "cog_name",
    ]
    assert df_mapping.to_dicts() == [
        {
            "name": "9606.ENSP1",
            "og": "OG0001",
            "level": "2759",
            "description": "alpha OG",
            "COG_categories": "EG",
            "cog_category": "E",
            "cog_class": "3",
            "cog_name": "Amino acid transport and metabolism",
        },
        {
            "name": "9606.ENSP1",
            "og": "OG0001",
            "level": "2759",
            "description": "alpha OG",
            "COG_categories": "EG",
            "cog_category": "G",
            "cog_class": "3",
            "cog_name": "Carbohydrate transport and metabolism",
        },
        {
            "name": "9606.ENSP1",
            "og": "OG0002",
            "level": "2759",
            "description": "beta OG",
            "COG_categories": "S",
            "cog_category": "S",
            "cog_class": "4",
            "cog_name": "Function unknown",
        },
    ]


def test_select_ids_queries_subset_and_reports_unmapped(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggNOGDatabase.from_sqlite(
        files["db"],
        cog_functions=files["cog_fun"],
    )

    selection = db.select_ids(
        ["9606.ENSP1", "9606.MISSING"],
    )

    df_mapping = selection.mappings().collect()
    assert df_mapping.select(
        "input_id", "input_namespace", "cog_category"
    ).to_dicts() == [
        {
            "input_id": "9606.ENSP1",
            "input_namespace": "eggnog_protein",
            "cog_category": "E",
        },
        {
            "input_id": "9606.ENSP1",
            "input_namespace": "eggnog_protein",
            "cog_category": "G",
        },
        {
            "input_id": "9606.ENSP1",
            "input_namespace": "eggnog_protein",
            "cog_category": "S",
        },
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"input_id": "9606.MISSING"}
    ]


def test_select_groups_preserves_group_id(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    db = EggNOGDatabase.from_sqlite(
        files["db"],
        cog_functions=files["cog_fun"],
    )

    selection = db.select_groups(
        {
            "up": ["9606.ENSP1"],
            "down": ["9606.ENSP1", "9606.MISSING"],
        },
    )

    df_mapping = selection.mappings().collect()
    assert selection.mappings().collect().equals(df_mapping)
    assert df_mapping.columns[:3] == [
        "group_id",
        "input_id",
        "input_namespace",
    ]
    assert df_mapping.select("group_id", "cog_category").to_dicts() == [
        {"group_id": "down", "cog_category": "E"},
        {"group_id": "down", "cog_category": "G"},
        {"group_id": "down", "cog_category": "S"},
        {"group_id": "up", "cog_category": "E"},
        {"group_id": "up", "cog_category": "G"},
        {"group_id": "up", "cog_category": "S"},
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "down", "input_id": "9606.MISSING"}
    ]


def test_grouped_selection_queries_one_global_unique_id_set(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    files = write_eggnog_fixture(tmp_path)
    queried_ids: list[list[str] | None] = []
    original = eggnog_util.iter_protein_ogs

    def counted_iter_protein_ogs(
        conn: sqlite3.Connection,
        eggnog_protein_ids: Iterable[str] | None,
    ) -> Iterator[tuple[str, str]]:
        queried = None if eggnog_protein_ids is None else list(eggnog_protein_ids)
        queried_ids.append(queried)
        return original(conn, queried)

    monkeypatch.setattr(eggnog_util, "iter_protein_ogs", counted_iter_protein_ogs)
    selection = EggNOGDatabase.from_sqlite(files["db"]).select_groups(
        {
            "up": ["9606.ENSP2", "9606.ENSP1", "9606.ENSP1"],
            "down": ["9606.ENSP1", "9606.MISSING"],
        }
    )

    mapping = selection.mappings().collect()
    selection.mappings().collect()
    selection.unmatched_ids().collect()
    selection.unmatched_ids().collect()

    assert mapping.height == 6
    expected_ids = ["9606.ENSP1", "9606.ENSP2", "9606.MISSING"]
    assert queried_ids
    assert all(ids == expected_ids for ids in queried_ids)


def test_gzip_sqlite_is_decompressed_to_tmp_dir(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)
    file_gz = tmp_path / "eggnog.db.gz"
    with (
        files["db"].open("rb") as handle_in,
        gzip.open(file_gz, "wb") as handle_out,
    ):
        shutil.copyfileobj(handle_in, handle_out)

    scratch = tmp_path / "tmp"
    with pytest.warns(UserWarning, match="Compressed eggNOG SQLite source detected"):
        db = EggNOGDatabase.from_sqlite(
            file_gz,
            cog_functions=files["cog_fun"],
            temp_dir=scratch,
        )

    assert (
        db.select_ids(
            ["9606.ENSP1"],
        )
        .mappings()
        .collect()
        .height
        == 3
    )
    assert list(scratch.iterdir()) == []
    assert {path.name for path in tmp_path.iterdir()} == {
        "cog-24.fun.tab",
        "eggnog.db",
        "eggnog.db.gz",
        "tmp",
    }


def test_gzip_sqlite_cleanup_on_query_failure(tmp_path: Path) -> None:
    file_gz = tmp_path / "broken.db.gz"
    with gzip.open(file_gz, "wb") as handle:
        handle.write(b"not a SQLite database")
    scratch = tmp_path / "scratch"

    with pytest.warns(UserWarning):
        db = EggNOGDatabase.from_sqlite(file_gz, temp_dir=scratch)
    with pytest.raises((sqlite3.DatabaseError, pl.exceptions.ComputeError)):
        db.select_ids(["9606.ENSP1"]).mappings().collect()

    assert list(scratch.iterdir()) == []
    assert {path.name for path in tmp_path.iterdir()} == {"broken.db.gz", "scratch"}


def test_plain_sqlite_emits_no_warning(tmp_path: Path) -> None:
    files = write_eggnog_fixture(tmp_path)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        db = EggNOGDatabase.from_sqlite(files["db"])
        db.select_ids(["9606.ENSP1"]).mappings().collect()
