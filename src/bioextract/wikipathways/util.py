from __future__ import annotations

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


def read_gmt_frames(file_gmt: Path) -> dict[str, pl.LazyFrame]:
    pathways: list[_PathwayRecord] = []
    term2gene_rows: list[_Term2GeneRecord] = []

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
            url = fields[1].strip()
            gene_ids = [gene_id.strip() for gene_id in fields[2:] if gene_id.strip()]
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
                        "WikiPathwaysId": pathway["WikiPathwaysId"],
                        "GeneId": gene_id,
                    }
                )

    lf_pathway = (
        pl.DataFrame(pathways, schema=SCHEMA_PATHWAY)
        .lazy()
        .unique(subset=["WikiPathwaysId"])
        .sort("WikiPathwaysId")
    )
    lf_term2gene = (
        pl.DataFrame(term2gene_rows, schema=SCHEMA_TERM2GENE)
        .lazy()
        .unique()
        .sort("WikiPathwaysId", "GeneId")
    )
    return {
        "pathway": lf_pathway,
        "term2gene": lf_term2gene,
        "term2name": extract_term2name_frame(lf_pathway),
    }


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
