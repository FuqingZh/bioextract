from __future__ import annotations

import glob as glob_module
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

import polars as pl

from .constant import SCHEMA_PATHWAY, SCHEMA_TERM2GENE, SCHEMA_TERM2NAME


class _PathwayRecord(TypedDict):
    PathwayName: str
    Collection: str
    Version: str
    WikiPathwaysId: str
    Species: str
    Url: str
    GeneCount: int


class _PathwayHeaderRecord(TypedDict):
    PathwayName: str
    Collection: str
    Version: str
    WikiPathwaysId: str
    Species: str


class _Term2GeneRecord(TypedDict):
    WikiPathwaysId: str
    GeneId: str


@dataclass(frozen=True, slots=True)
class _ParsedGmtDataset:
    frames: dict[str, pl.LazyFrame]
    release_version: str


def resolve_gmt_sources(
    source: os.PathLike[str] | str | Sequence[os.PathLike[str] | str],
    *,
    glob: bool,
) -> tuple[Path, ...]:
    entries: tuple[os.PathLike[str] | str, ...] = (
        (source,) if isinstance(source, (str, os.PathLike)) else tuple(source)
    )
    if not entries:
        raise ValueError("WikiPathways GMT source must contain at least one path")

    matched: list[Path] = []
    for entry in entries:
        try:
            value = os.fspath(entry)
        except TypeError as error:
            raise TypeError(
                "WikiPathways GMT source entries must be strings or path-like"
            ) from error
        if not value:
            raise ValueError("WikiPathways GMT source paths must be non-empty")
        candidates = (
            [Path(match) for match in glob_module.glob(value, recursive=True)]
            if glob
            else [Path(value)]
        )
        if not candidates:
            raise FileNotFoundError(
                f"WikiPathways GMT source pattern matched no files: {value}"
            )
        matched.extend(candidates)

    resolved: list[Path] = []
    seen: set[Path] = set()
    for candidate in matched:
        try:
            actual = candidate.resolve(strict=True)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"WikiPathways GMT file not found: {candidate}"
            ) from error
        if not actual.is_file():
            raise ValueError(f"WikiPathways GMT source is not a file: {candidate}")
        if actual in seen:
            raise ValueError(
                "WikiPathways GMT source resolves to a duplicate physical file: "
                f"{actual}"
            )
        seen.add(actual)
        resolved.append(actual)
    return tuple(sorted(resolved, key=os.fspath))


def read_gmt_frames(files_gmt: tuple[Path, ...]) -> _ParsedGmtDataset:
    pathways: list[_PathwayRecord] = []
    term2gene_rows: list[_Term2GeneRecord] = []

    pathway_locations: dict[str, tuple[Path, int]] = {}
    for file_gmt in files_gmt:
        with file_gmt.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"WikiPathways GMT line must contain at least two columns: "
                        f"path={file_gmt}, line={line_number}"
                    )
                pathway = parse_pathway_header(
                    fields[0], file_gmt=file_gmt, line_number=line_number
                )
                pathway_id = pathway["WikiPathwaysId"]
                if pathway_id in pathway_locations:
                    first_path, first_line = pathway_locations[pathway_id]
                    raise ValueError(
                        "WikiPathwaysId must be unique across all GMT files: "
                        f"id={pathway_id}, first_path={first_path}, "
                        f"first_line={first_line}, duplicate_path={file_gmt}, "
                        f"duplicate_line={line_number}"
                    )
                pathway_locations[pathway_id] = (file_gmt, line_number)
                url = fields[1].strip()
                gene_ids = [
                    gene_id.strip() for gene_id in fields[2:] if gene_id.strip()
                ]
                gene_ids_unique = sorted(set(gene_ids))
                pathway_record: _PathwayRecord = {
                    **pathway,
                    "Url": url,
                    "GeneCount": len(gene_ids_unique),
                }
                pathways.append(pathway_record)
                for gene_id in gene_ids_unique:
                    term2gene_rows.append(
                        {
                            "WikiPathwaysId": pathway_id,
                            "GeneId": gene_id,
                        }
                    )

    collections = {pathway["Collection"] for pathway in pathways}
    versions = {pathway["Version"] for pathway in pathways}
    if len(collections) != 1:
        raise ValueError(
            "WikiPathways GMT files must contain one common Collection; "
            f"found={sorted(collections)}"
        )
    if len(versions) != 1:
        raise ValueError(
            "WikiPathways GMT files must contain one common Version; "
            f"found={sorted(versions)}"
        )

    lf_pathway = (
        pl.DataFrame(pathways, schema=SCHEMA_PATHWAY).lazy().sort("WikiPathwaysId")
    )
    lf_term2gene = (
        pl.DataFrame(term2gene_rows, schema=SCHEMA_TERM2GENE)
        .lazy()
        .unique()
        .sort("WikiPathwaysId", "GeneId")
    )
    return _ParsedGmtDataset(
        frames={
            "pathway": lf_pathway,
            "term2gene": lf_term2gene,
            "term2name": extract_term2name_frame(lf_pathway),
        },
        release_version=next(iter(versions)),
    )


def parse_pathway_header(
    header: str,
    *,
    file_gmt: Path,
    line_number: int,
) -> _PathwayHeaderRecord:
    parts = [part.strip() for part in header.split("%")]
    if len(parts) != 4 or any(not part for part in parts):
        raise ValueError(
            "WikiPathways GMT header must have four '%' separated fields: "
            "PathwayName%Collection%WikiPathwaysId%Species; "
            f"path={file_gmt}, line={line_number}, value={header!r}"
        )
    pathway_name, collection, pathway_id, species = parts
    version = collection.removeprefix("WikiPathways_")
    return {
        "PathwayName": pathway_name,
        "Collection": collection,
        "Version": version,
        "WikiPathwaysId": pathway_id,
        "Species": species,
    }


def extract_term2name_frame(lf_pathway: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf_pathway.select(SCHEMA_TERM2NAME.keys())
        .unique(subset=["WikiPathwaysId"])
        .sort("WikiPathwaysId")
    )


def extract_mapping_frame(
    lf_pathway: pl.LazyFrame,
    lf_term2gene: pl.LazyFrame,
    df_input_ids: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.LazyFrame:
    cols_group = list(cols_group_id)
    cols_out = cols_group + [
        "InputId",
        "GeneId",
        "WikiPathwaysId",
        "PathwayName",
        "Species",
        "Url",
    ]
    return (
        df_input_ids.lazy()
        .join(
            lf_term2gene,
            left_on="InputId",
            right_on="GeneId",
            how="inner",
        )
        .with_columns(pl.col("InputId").alias("GeneId"))
        .join(
            lf_pathway.select("WikiPathwaysId", "PathwayName", "Species", "Url"),
            on="WikiPathwaysId",
            how="left",
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )


def extract_unmatched_ids_frame(
    df_input_ids: pl.DataFrame,
    lf_mapping: pl.LazyFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.LazyFrame:
    cols_index = list(cols_group_id) + ["InputId"]
    lf_mapped_input_ids = lf_mapping.select(cols_index).unique().sort(cols_index)
    return (
        df_input_ids.lazy()
        .join(lf_mapped_input_ids, on=cols_index, how="anti")
        .select(cols_index)
        .sort(cols_index)
    )
