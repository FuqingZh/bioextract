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
    pathway_name: str
    collection: str
    version: str
    wiki_pathways_id: str
    species: str
    url: str
    gene_count: int


class _PathwayHeaderRecord(TypedDict):
    pathway_name: str
    collection: str
    version: str
    wiki_pathways_id: str
    species: str


class _Term2GeneRecord(TypedDict):
    wiki_pathways_id: str
    gene_id: str


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
    seen_paths: set[Path] = set()
    seen_file_ids: set[tuple[int, int]] = set()
    for candidate in matched:
        try:
            actual = candidate.resolve(strict=True)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"WikiPathways GMT file not found: {candidate}"
            ) from error
        if not actual.is_file():
            raise ValueError(f"WikiPathways GMT source is not a file: {candidate}")
        stat = actual.stat()
        file_id = (stat.st_dev, stat.st_ino)
        if actual in seen_paths or file_id in seen_file_ids:
            raise ValueError(
                "WikiPathways GMT source resolves to a duplicate physical file: "
                f"{actual}"
            )
        seen_paths.add(actual)
        seen_file_ids.add(file_id)
        resolved.append(actual)
    return tuple(sorted(resolved, key=os.fspath))


def read_gmt_frames(files_gmt: tuple[Path, ...]) -> _ParsedGmtDataset:
    pathways: list[_PathwayRecord] = []
    term2gene_rows: list[_Term2GeneRecord] = []

    pathway_locations: dict[str, tuple[Path, int]] = {}
    for file_gmt in files_gmt:
        pathway_count = 0
        with file_gmt.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.rstrip("\n")
                if not line.strip():
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
                pathway_id = pathway["wiki_pathways_id"]
                if pathway_id in pathway_locations:
                    first_path, first_line = pathway_locations[pathway_id]
                    raise ValueError(
                        "wiki_pathways_id must be unique across all GMT files: "
                        f"id={pathway_id}, first_path={first_path}, "
                        f"first_line={first_line}, duplicate_path={file_gmt}, "
                        f"duplicate_line={line_number}"
                    )
                pathway_locations[pathway_id] = (file_gmt, line_number)
                pathway_count += 1
                url = fields[1].strip()
                gene_ids = [
                    gene_id.strip() for gene_id in fields[2:] if gene_id.strip()
                ]
                gene_ids_unique = sorted(set(gene_ids))
                pathway_record: _PathwayRecord = {
                    **pathway,
                    "url": url,
                    "gene_count": len(gene_ids_unique),
                }
                pathways.append(pathway_record)
                for gene_id in gene_ids_unique:
                    term2gene_rows.append(
                        {
                            "wiki_pathways_id": pathway_id,
                            "gene_id": gene_id,
                        }
                    )
        if pathway_count == 0:
            raise ValueError(
                "WikiPathways GMT file must contain at least one non-empty "
                f"pathway record: path={file_gmt}"
            )

    collections = {pathway["collection"] for pathway in pathways}
    versions = {pathway["version"] for pathway in pathways}
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
        pl.DataFrame(pathways, schema=SCHEMA_PATHWAY).lazy().sort("wiki_pathways_id")
    )
    lf_term2gene = (
        pl.DataFrame(term2gene_rows, schema=SCHEMA_TERM2GENE)
        .lazy()
        .unique()
        .sort("wiki_pathways_id", "gene_id")
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
    collection_prefix = "WikiPathways_"
    if not collection.startswith(collection_prefix):
        raise ValueError(
            "WikiPathways GMT Collection must start with 'WikiPathways_': "
            f"path={file_gmt}, line={line_number}, value={collection!r}"
        )
    version = collection[len(collection_prefix) :]
    if not version:
        raise ValueError(
            "WikiPathways GMT Collection must include a non-empty Version after "
            f"'WikiPathways_': path={file_gmt}, line={line_number}"
        )
    return {
        "pathway_name": pathway_name,
        "collection": collection,
        "version": version,
        "wiki_pathways_id": pathway_id,
        "species": species,
    }


def extract_term2name_frame(lf_pathway: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf_pathway.select(SCHEMA_TERM2NAME.keys())
        .unique(subset=["wiki_pathways_id"])
        .sort("wiki_pathways_id")
    )


def extract_mapping_frame(
    lf_pathway: pl.LazyFrame,
    lf_term2gene: pl.LazyFrame,
    df_input_ids: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.LazyFrame:
    cols_out = [
        "input_id",
        "gene_id",
        "wiki_pathways_id",
        "pathway_name",
        "species",
        "url",
    ]
    lf_hits = (
        df_input_ids.lazy()
        .join(
            lf_term2gene,
            left_on="input_id",
            right_on="gene_id",
            how="inner",
        )
        .with_columns(pl.col("input_id").alias("gene_id"))
        .join(
            lf_pathway.select("wiki_pathways_id", "pathway_name", "species", "url"),
            on="wiki_pathways_id",
            how="left",
        )
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )
    if df_group_membership is None:
        return lf_hits
    grouped_cols = ["group_id", *cols_out]
    return (
        df_group_membership.lazy()
        .join(lf_hits, on="input_id", how="inner")
        .select(grouped_cols)
        .unique()
        .sort(grouped_cols)
    )


def extract_unmatched_ids_frame(
    df_input_ids: pl.DataFrame,
    df_mapping: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    df_mapped_input_ids = df_mapping.select("input_id").unique().sort("input_id")
    df_unmatched = (
        df_input_ids.join(df_mapped_input_ids, on="input_id", how="anti")
        .select("input_id")
        .sort("input_id")
    )
    if df_group_membership is None:
        return df_unmatched
    return (
        df_group_membership.join(df_unmatched, on="input_id", how="inner")
        .select("group_id", "input_id")
        .sort("group_id", "input_id")
    )
