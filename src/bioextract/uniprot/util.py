from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO

import polars as pl

from bioextract._shared import validate_required_cols

from .constant import (
    COLS_EGGNOG_XREF,
    COLS_IDMAPPING_SELECTED,
    REQUIRED_COLS_MAPPING,
    SCHEMA_EGGNOG_XREF,
    SCHEMA_MAPPING,
)


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


def has_hive_parquet_candidates(dir_mapping: Path) -> bool:
    return any(path.is_file() for path in dir_mapping.glob("**/*.parquet"))


def read_eggnog_xref_frame(
    file_dat: Path,
    *,
    source_db: str,
    input_ids: set[str] | None = None,
) -> pl.DataFrame:
    rows: list[dict[str, str | bool]] = []
    with open_uniprot_dat(file_dat) as handle:
        accessions: list[str] = []
        eggnog_xrefs: list[tuple[str, str]] = []
        for line in handle:
            if line.startswith("//"):
                append_eggnog_xref_rows(
                    rows,
                    accessions=accessions,
                    eggnog_xrefs=eggnog_xrefs,
                    source_db=source_db,
                    input_ids=input_ids,
                )
                accessions = []
                eggnog_xrefs = []
                continue
            if line.startswith("AC   "):
                accessions.extend(parse_accession_line(line))
            elif line.startswith("DR   eggNOG;"):
                eggnog_xrefs.append(parse_eggnog_dr_line(line))

        append_eggnog_xref_rows(
            rows,
            accessions=accessions,
            eggnog_xrefs=eggnog_xrefs,
            source_db=source_db,
            input_ids=input_ids,
        )

    if not rows:
        return pl.DataFrame(schema=SCHEMA_EGGNOG_XREF)
    return pl.DataFrame(rows, schema=SCHEMA_EGGNOG_XREF).unique().sort(COLS_EGGNOG_XREF)


def write_eggnog_xref_tsv(
    file_dat: Path,
    file_out: Path,
    *,
    source_db: str,
) -> None:
    file_out.parent.mkdir(parents=True, exist_ok=True)
    with file_out.open("w", encoding="utf-8") as handle_out:
        handle_out.write("\t".join(COLS_EGGNOG_XREF) + "\n")
        with open_uniprot_dat(file_dat) as handle_in:
            accessions: list[str] = []
            eggnog_xrefs: list[tuple[str, str]] = []
            for line in handle_in:
                if line.startswith("//"):
                    write_eggnog_xref_rows(
                        handle_out,
                        accessions=accessions,
                        eggnog_xrefs=eggnog_xrefs,
                        source_db=source_db,
                    )
                    accessions = []
                    eggnog_xrefs = []
                    continue
                if line.startswith("AC   "):
                    accessions.extend(parse_accession_line(line))
                elif line.startswith("DR   eggNOG;"):
                    eggnog_xrefs.append(parse_eggnog_dr_line(line))

            write_eggnog_xref_rows(
                handle_out,
                accessions=accessions,
                eggnog_xrefs=eggnog_xrefs,
                source_db=source_db,
            )


def scan_eggnog_xref_tsv(file_tsv: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_tsv,
        separator="\t",
        has_header=True,
        schema_overrides=SCHEMA_EGGNOG_XREF,
    ).select(COLS_EGGNOG_XREF)


def append_eggnog_xref_rows(
    rows: list[dict[str, str | bool]],
    *,
    accessions: list[str],
    eggnog_xrefs: list[tuple[str, str]],
    source_db: str,
    input_ids: set[str] | None,
) -> None:
    if not accessions or not eggnog_xrefs:
        return
    primary_accession = accessions[0]
    for accession in accessions:
        if input_ids is not None and accession not in input_ids:
            continue
        for eggnog_og_id, eggnog_level in eggnog_xrefs:
            rows.append(
                {
                    "UniProtId": accession,
                    "PrimaryUniProtId": primary_accession,
                    "IsPrimaryAccession": accession == primary_accession,
                    "EggnogOgId": eggnog_og_id,
                    "EggnogLevel": eggnog_level,
                    "SourceDb": source_db,
                }
            )


def write_eggnog_xref_rows(
    handle: TextIO,
    *,
    accessions: list[str],
    eggnog_xrefs: list[tuple[str, str]],
    source_db: str,
) -> None:
    if not accessions or not eggnog_xrefs:
        return
    primary_accession = accessions[0]
    for accession in accessions:
        is_primary = "true" if accession == primary_accession else "false"
        for eggnog_og_id, eggnog_level in eggnog_xrefs:
            handle.write(
                "\t".join(
                    [
                        accession,
                        primary_accession,
                        is_primary,
                        eggnog_og_id,
                        eggnog_level,
                        source_db,
                    ]
                )
                + "\n"
            )


def parse_accession_line(line: str) -> list[str]:
    return [
        accession.strip()
        for accession in line[5:].strip().split(";")
        if accession.strip()
    ]


def parse_eggnog_dr_line(line: str) -> tuple[str, str]:
    payload = line.removeprefix("DR   eggNOG;").strip()
    parts = [part.strip().rstrip(".") for part in payload.split(";") if part.strip()]
    if len(parts) != 2:
        raise ValueError(f"Unsupported UniProt eggNOG DR line: {line.rstrip()!r}")
    return parts[0], parts[1]


def open_uniprot_dat(file_dat: Path) -> TextIO:
    if file_dat.name.endswith(".gz"):
        return gzip.open(file_dat, "rt", encoding="utf-8", errors="replace")
    return file_dat.open("r", encoding="utf-8", errors="replace")
