from __future__ import annotations

from pathlib import Path

import polars as pl

from bioextract._shared import validate_required_cols

from .constant import COLS_IDMAPPING_SELECTED, REQUIRED_COLS_MAPPING, SCHEMA_MAPPING


def normalize_taxids(taxids: tuple[str | int, ...]) -> tuple[str, ...]:
    taxids_normalized = tuple(
        dict.fromkeys(str(taxid).strip() for taxid in taxids if str(taxid).strip())
    )
    if len(taxids_normalized) != len(taxids):
        raise ValueError("TaxId values must be non-empty after normalization")
    return taxids_normalized


def scan_raw_idmapping_selected(file_idmapping_selected: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_idmapping_selected,
        separator="\t",
        has_header=False,
        new_columns=COLS_IDMAPPING_SELECTED,
        schema_overrides=SCHEMA_MAPPING,
        infer_schema=False,
        quote_char=None,
        missing_utf8_is_empty_string=True,
    ).select(COLS_IDMAPPING_SELECTED)


def scan_parquet_mapping(file_mapping: Path) -> pl.LazyFrame:
    return pl.scan_parquet(file_mapping)


def scan_hive_mapping_dataset(dir_mapping: Path) -> pl.LazyFrame:
    return pl.scan_parquet(
        dir_mapping / "**" / "*.parquet",
        hive_partitioning=True,
        hive_schema={"TaxId": pl.String},
    )


def filter_taxids(lf: pl.LazyFrame, taxids: tuple[str, ...]) -> pl.LazyFrame:
    if not taxids:
        return lf
    return lf.filter(pl.col("TaxId").is_in(taxids))


def validate_mapping_schema(lf: pl.LazyFrame) -> None:
    schema = lf.collect_schema()
    validate_required_cols(
        cols_available=schema.names(),
        cols_required=REQUIRED_COLS_MAPPING,
        context="UniProt idmapping selected mapping",
    )


def has_parquet_files(dir_mapping: Path) -> bool:
    return any(
        path.is_file() and path.suffix == ".parquet" for path in dir_mapping.rglob("*")
    )
