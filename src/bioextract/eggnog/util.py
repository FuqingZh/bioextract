from __future__ import annotations

import csv
import gzip
import shutil
import sqlite3
import tempfile
from collections.abc import Generator, Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Protocol
from typing import Sequence

import polars as pl

from .constant import (
    COLS_MAPPING,
    KIND_INPUT_ID_VALUES,
    SCHEMA_MAPPING,
    EggnogInputIdKind,
)


class _RowWriter(Protocol):
    def writerow(self, row: Sequence[object]) -> object: ...


def read_cog_fun_frame(file_cog_fun: Path | None) -> pl.DataFrame:
    schema = {
        "CogCategory": pl.String,
        "CogClass": pl.String,
        "_Color": pl.String,
        "CogName": pl.String,
    }
    if file_cog_fun is None:
        return pl.DataFrame(
            schema={
                "CogCategory": pl.String,
                "CogClass": pl.String,
                "CogName": pl.String,
            }
        )
    if file_cog_fun.stat().st_size == 0:
        return pl.DataFrame(schema=schema).select("CogCategory", "CogClass", "CogName")
    df = pl.scan_csv(
        file_cog_fun,
        separator="\t",
        has_header=False,
        new_columns=["CogCategory", "CogClass", "_Color", "CogName"],
        schema_overrides=schema,
        truncate_ragged_lines=True,
    ).collect()
    return (
        df.filter(pl.col("CogCategory").str.len_chars() == 1)
        .select("CogCategory", "CogClass", "CogName")
        .unique()
        .sort("CogCategory")
    )


def scan_mapping_tsv(file_mapping_tsv: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_mapping_tsv,
        separator="\t",
        has_header=True,
        schema_overrides=SCHEMA_MAPPING,
    ).select(COLS_MAPPING)


def write_mapping_tsv(
    *,
    file_eggnog_db: Path,
    dir_tmp: Path | None,
    df_cog_fun: pl.DataFrame,
    file_out: Path,
) -> None:
    file_out.parent.mkdir(parents=True, exist_ok=True)
    map_cog_fun = {
        row["CogCategory"]: (row["CogClass"], row["CogName"])
        for row in df_cog_fun.iter_rows(named=True)
    }
    with open_sqlite_path(file_eggnog_db, dir_tmp=dir_tmp) as file_sqlite:
        with sqlite3.connect(file_sqlite) as conn:
            with file_out.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(
                    handle,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writerow(COLS_MAPPING)
                write_mapping_rows(writer, conn=conn, map_cog_fun=map_cog_fun)


def write_mapping_rows(
    writer: _RowWriter,
    *,
    conn: sqlite3.Connection,
    map_cog_fun: dict[str, tuple[str | None, str | None]],
) -> None:
    og_cache: dict[tuple[str, str | None], list[dict[str, str | None]]] = {}
    for protein_id, ogs_text in iter_protein_ogs(conn, None):
        for og_id, og_level in parse_ogs(ogs_text):
            key = (og_id, og_level)
            if key not in og_cache:
                og_cache[key] = read_og_rows(conn, [key])
            for og_row in og_cache[key]:
                for category in parse_cog_categories(og_row["CogCategories"]):
                    cog_class, cog_name = map_cog_fun.get(category, (None, None))
                    writer.writerow(
                        [
                            protein_id,
                            og_row["EggnogOgId"] or "",
                            og_row["EggnogLevel"] or "",
                            category,
                            cog_class or "",
                            cog_name or "",
                            og_row["OgDescription"] or "",
                        ]
                    )


def build_mapping_frame(
    *,
    file_eggnog_db: Path,
    dir_tmp: Path | None,
    df_cog_fun: pl.DataFrame,
) -> pl.DataFrame:
    with open_sqlite_path(file_eggnog_db, dir_tmp=dir_tmp) as file_sqlite:
        return read_mapping_frame_from_sqlite(file_sqlite, df_cog_fun=df_cog_fun)


def select_mapping_frame(
    *,
    file_eggnog_db: Path,
    dir_tmp: Path | None,
    df_input_ids: pl.DataFrame,
    kind_input_id: EggnogInputIdKind,
    cols_group_id: tuple[str, ...],
    df_cog_fun: pl.DataFrame,
) -> pl.DataFrame:
    validate_kind_input_id(kind_input_id)
    if df_input_ids.height == 0:
        cols_out = list(cols_group_id) + ["InputId", "KindInputId"] + COLS_MAPPING
        return pl.DataFrame(schema={col: pl.String for col in cols_out})

    input_ids = df_input_ids.get_column("InputId").unique().sort().to_list()
    with open_sqlite_path(file_eggnog_db, dir_tmp=dir_tmp) as file_sqlite:
        df_mapping = read_mapping_frame_from_sqlite(
            file_sqlite,
            df_cog_fun=df_cog_fun,
            eggnog_protein_ids=input_ids,
        )

    return extract_mapping_frame(
        df_mapping,
        df_input_ids,
        kind_input_id=kind_input_id,
        cols_group_id=cols_group_id,
    )


def read_mapping_frame_from_sqlite(
    file_sqlite: Path,
    *,
    df_cog_fun: pl.DataFrame,
    eggnog_protein_ids: Iterable[str] | None = None,
) -> pl.DataFrame:
    rows: list[dict[str, str | None]] = []
    with sqlite3.connect(file_sqlite) as conn:
        for protein_id, ogs_text in iter_protein_ogs(conn, eggnog_protein_ids):
            for og_id, og_level in parse_ogs(ogs_text):
                rows.append(
                    {
                        "EggnogProteinId": protein_id,
                        "EggnogOgId": og_id,
                        "EggnogLevel": og_level,
                    }
                )

        if not rows:
            return pl.DataFrame(schema=SCHEMA_MAPPING)

        df_protein_og = pl.DataFrame(
            rows,
            schema={
                "EggnogProteinId": pl.String,
                "EggnogOgId": pl.String,
                "EggnogLevel": pl.String,
            },
        ).unique()

        og_rows = read_og_rows(
            conn,
            df_protein_og.select("EggnogOgId", "EggnogLevel").unique().iter_rows(),
        )

    if not og_rows:
        return pl.DataFrame(schema=SCHEMA_MAPPING)

    df_og = pl.DataFrame(
        og_rows,
        schema={
            "EggnogOgId": pl.String,
            "EggnogLevel": pl.String,
            "OgDescription": pl.String,
            "CogCategories": pl.String,
        },
    )
    df_categories = expand_cog_categories(df_og)
    if df_categories.height == 0:
        return pl.DataFrame(schema=SCHEMA_MAPPING)

    return (
        df_protein_og.join(
            df_categories,
            on=["EggnogOgId", "EggnogLevel"],
            how="inner",
        )
        .join(df_cog_fun, on="CogCategory", how="left")
        .select(COLS_MAPPING)
        .unique()
        .sort(COLS_MAPPING)
    )


def iter_protein_ogs(
    conn: sqlite3.Connection,
    eggnog_protein_ids: Iterable[str] | None,
) -> Iterator[tuple[str, str]]:
    if eggnog_protein_ids is None:
        cursor = conn.execute("select name, ogs from prots where ogs is not null")
        yield from ((str(name), str(ogs)) for name, ogs in cursor if str(ogs))
        return

    protein_ids = [str(protein_id) for protein_id in eggnog_protein_ids]
    for chunk in chunked(protein_ids, 500):
        placeholders = ",".join("?" for _ in chunk)
        cursor = conn.execute(
            f"select name, ogs from prots where name in ({placeholders}) "
            "and ogs is not null",
            chunk,
        )
        yield from ((str(name), str(ogs)) for name, ogs in cursor if str(ogs))


def read_og_rows(
    conn: sqlite3.Connection,
    og_keys: Iterable[tuple[str, str | None]],
) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    keys_by_level: dict[str, list[str]] = {}
    keys_without_level: list[str] = []
    for og_id, og_level in og_keys:
        if og_level:
            keys_by_level.setdefault(og_level, []).append(og_id)
        else:
            keys_without_level.append(og_id)

    for og_level, og_ids in keys_by_level.items():
        for chunk in chunked(sorted(set(og_ids)), 500):
            placeholders = ",".join("?" for _ in chunk)
            cursor = conn.execute(
                "select og, level, description, COG_categories from og "
                f"where level = ? and og in ({placeholders})",
                [og_level, *chunk],
            )
            rows.extend(format_og_query_rows(cursor))

    for chunk in chunked(sorted(set(keys_without_level)), 500):
        placeholders = ",".join("?" for _ in chunk)
        cursor = conn.execute(
            "select og, level, description, COG_categories from og "
            f"where og in ({placeholders})",
            chunk,
        )
        rows.extend(format_og_query_rows(cursor))

    return rows


def format_og_query_rows(
    cursor: sqlite3.Cursor,
) -> list[dict[str, str | None]]:
    return [
        {
            "EggnogOgId": str(og_id),
            "EggnogLevel": str(level) if level is not None else None,
            "OgDescription": str(description) if description is not None else None,
            "CogCategories": str(categories) if categories is not None else None,
        }
        for og_id, level, description, categories in cursor
    ]


def parse_ogs(value: str | None) -> list[tuple[str, str | None]]:
    if value is None:
        return []
    pairs: list[tuple[str, str | None]] = []
    for item in str(value).replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        if "@" in item:
            og_id, level = item.split("@", 1)
            pairs.append((og_id.strip(), level.strip() or None))
        else:
            pairs.append((item, None))
    return pairs


def expand_cog_categories(df_og: pl.DataFrame) -> pl.DataFrame:
    rows: list[dict[str, str | None]] = []
    for row in df_og.iter_rows(named=True):
        for category in parse_cog_categories(row["CogCategories"]):
            rows.append(
                {
                    "EggnogOgId": row["EggnogOgId"],
                    "EggnogLevel": row["EggnogLevel"],
                    "OgDescription": row["OgDescription"],
                    "CogCategory": category,
                }
            )
    if not rows:
        return pl.DataFrame(
            schema={
                "EggnogOgId": pl.String,
                "EggnogLevel": pl.String,
                "OgDescription": pl.String,
                "CogCategory": pl.String,
            }
        )
    return pl.DataFrame(rows).unique()


def parse_cog_categories(value: str | None) -> list[str]:
    if value is None:
        return []
    return [
        category
        for category in dict.fromkeys(str(value).strip())
        if category and category != "-"
    ]


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    kind_input_id: EggnogInputIdKind,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    validate_kind_input_id(kind_input_id)
    cols_group = list(cols_group_id)
    cols_out = cols_group + ["InputId", "KindInputId"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="InputId",
            right_on="EggnogProteinId",
            how="inner",
        )
        .with_columns(
            pl.col("InputId").alias("EggnogProteinId"),
            pl.lit(kind_input_id).alias("KindInputId"),
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def extract_unmapped_input_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
    cols_index = list(cols_group_id) + ["InputId"]
    df_mapped_input_ids = df_mapping.select(cols_index).unique().sort(cols_index)
    return (
        df_input_ids.join(df_mapped_input_ids, on=cols_index, how="anti")
        .select(cols_index)
        .sort(cols_index)
    )


def validate_kind_input_id(kind_input_id: str) -> None:
    if kind_input_id not in KIND_INPUT_ID_VALUES:
        raise ValueError(
            "kind_input_id must be one of: "
            f"{', '.join(KIND_INPUT_ID_VALUES)}; got {kind_input_id!r}"
        )


@contextmanager
def open_sqlite_path(
    file_eggnog_db: Path,
    *,
    dir_tmp: Path | None,
) -> Generator[Path, None, None]:
    if file_eggnog_db.suffix != ".gz":
        yield file_eggnog_db
        return

    with tempfile.TemporaryDirectory(dir=dir_tmp) as dir_work:
        file_sqlite = Path(dir_work) / file_eggnog_db.with_suffix("").name
        with gzip.open(file_eggnog_db, "rb") as handle_in:
            with file_sqlite.open("wb") as handle_out:
                shutil.copyfileobj(handle_in, handle_out, length=1024 * 1024 * 16)
        yield file_sqlite


def chunked(values: list[str], size: int) -> Iterator[list[str]]:
    for index in range(0, len(values), size):
        yield values[index : index + size]
