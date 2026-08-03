from __future__ import annotations

import gzip
import shutil
import sqlite3
import tempfile
from collections.abc import Generator, Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path

import polars as pl

from bioextract._shared import RowWriter, create_tsv_writer

from .constant import (
    COLS_MAPPING,
    NAMESPACE_VALUES,
    SCHEMA_MAPPING,
    EggnogNamespace,
)


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
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    map_cog_fun = {
        row["CogCategory"]: (row["CogClass"], row["CogName"])
        for row in df_cog_fun.iter_rows(named=True)
    }
    with (
        open_sqlite_path(file_eggnog_db, dir_tmp=dir_tmp) as file_sqlite,
        sqlite3.connect(file_sqlite) as conn,
        path.open("w", encoding="utf-8", newline="") as handle,
    ):
        writer = create_tsv_writer(handle)
        writer.writerow(COLS_MAPPING)
        write_mapping_rows(writer, conn=conn, map_cog_fun=map_cog_fun)


def write_mapping_rows(
    writer: RowWriter,
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
    df_group_membership: pl.DataFrame | None,
    namespace: EggnogNamespace,
    df_cog_fun: pl.DataFrame,
) -> pl.DataFrame:
    validate_namespace(namespace)
    if df_input_ids.height == 0:
        cols_group = ["GroupId"] if df_group_membership is not None else []
        cols_out = cols_group + ["InputId", "InputNamespace"] + COLS_MAPPING
        return pl.DataFrame(schema=dict.fromkeys(cols_out, pl.String))

    input_ids = df_input_ids.get_column("InputId").unique().sort().to_list()
    with open_sqlite_path(file_eggnog_db, dir_tmp=dir_tmp) as file_sqlite:
        df_mapping = read_mapping_frame_from_sqlite(
            file_sqlite,
            df_cog_fun=df_cog_fun,
            eggnog_protein_ids=input_ids,
        )

    df_selected = extract_mapping_frame(
        df_mapping,
        df_input_ids,
        namespace=namespace,
    )
    return _expand_group_membership(df_selected, df_group_membership)


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
    namespace: EggnogNamespace,
) -> pl.DataFrame:
    validate_namespace(namespace)
    cols_out = ["InputId", "InputNamespace"] + COLS_MAPPING
    return (
        df_input_ids.join(
            df_mapping,
            left_on="InputId",
            right_on="EggnogProteinId",
            how="inner",
        )
        .with_columns(
            pl.col("InputId").alias("EggnogProteinId"),
            pl.lit(namespace).alias("InputNamespace"),
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def extract_unmatched_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    df_mapped_input_ids = df_mapping.select("InputId").unique().sort("InputId")
    df_unmatched = (
        df_input_ids.join(df_mapped_input_ids, on="InputId", how="anti")
        .select("InputId")
        .sort("InputId")
    )
    if df_group_membership is None:
        return df_unmatched
    return (
        df_group_membership.join(df_unmatched, on="InputId", how="inner")
        .select("GroupId", "InputId")
        .sort("GroupId", "InputId")
    )


def _expand_group_membership(
    df_mapping: pl.DataFrame,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    if df_group_membership is None:
        return df_mapping
    cols_out = ["GroupId", *df_mapping.columns]
    return (
        df_group_membership.join(df_mapping, on="InputId", how="inner")
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def validate_namespace(namespace: str) -> None:
    if namespace not in NAMESPACE_VALUES:
        raise ValueError(
            "namespace must be one of: "
            f"{', '.join(NAMESPACE_VALUES)}; got {namespace!r}"
        )


@contextmanager
def open_sqlite_path(
    file_eggnog_db: Path,
    *,
    dir_tmp: Path | None,
) -> Generator[Path]:
    if file_eggnog_db.suffix != ".gz":
        yield file_eggnog_db
        return

    with tempfile.TemporaryDirectory(dir=dir_tmp) as dir_work:
        file_sqlite = Path(dir_work) / file_eggnog_db.with_suffix("").name
        with (
            gzip.open(file_eggnog_db, "rb") as handle_in,
            file_sqlite.open("wb") as handle_out,
        ):
            shutil.copyfileobj(handle_in, handle_out, length=1024 * 1024 * 16)
        yield file_sqlite


def chunked(values: list[str], size: int) -> Iterator[list[str]]:
    for index in range(0, len(values), size):
        yield values[index : index + size]
