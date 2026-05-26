from __future__ import annotations

from pathlib import Path

import polars as pl

from .constant import SCHEMA_PATHWAY, SCHEMA_TERM2GENE, SCHEMA_TERM2NAME


def read_gmt_frames(file_gmt: Path) -> dict[str, pl.DataFrame]:
    pathways: list[dict[str, object]] = []
    term2gene_rows: list[dict[str, str]] = []

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
            pathway = parse_pathway_header(fields[0], file_gmt=file_gmt, line_number=line_number)
            url = fields[1].strip()
            gene_ids = [gene_id.strip() for gene_id in fields[2:] if gene_id.strip()]
            gene_ids_unique = sorted(set(gene_ids))
            pathways.append(
                {
                    **pathway,
                    "Url": url,
                    "GeneCount": len(gene_ids_unique),
                }
            )
            for gene_id in gene_ids_unique:
                term2gene_rows.append(
                    {
                        "WikiPathwaysId": pathway["WikiPathwaysId"],
                        "GeneId": gene_id,
                    }
                )

    df_pathway = pl.DataFrame(pathways, schema=SCHEMA_PATHWAY)
    df_term2gene = pl.DataFrame(term2gene_rows, schema=SCHEMA_TERM2GENE)
    if df_pathway.height > 0:
        df_pathway = df_pathway.unique(subset=["WikiPathwaysId"]).sort("WikiPathwaysId")
    if df_term2gene.height > 0:
        df_term2gene = df_term2gene.unique().sort("WikiPathwaysId", "GeneId")
    return {
        "pathway": df_pathway,
        "term2gene": df_term2gene,
        "term2name": extract_term2name_frame(df_pathway),
    }


def parse_pathway_header(
    header: str,
    *,
    file_gmt: Path,
    line_number: int,
) -> dict[str, str]:
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


def extract_term2name_frame(df_pathway: pl.DataFrame) -> pl.DataFrame:
    return (
        df_pathway.select(SCHEMA_TERM2NAME.keys())
        .unique(subset=["WikiPathwaysId"])
        .sort("WikiPathwaysId")
    )


def extract_mapping_frame(
    df_pathway: pl.DataFrame,
    df_term2gene: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    cols_group_id: tuple[str, ...],
) -> pl.DataFrame:
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
        df_input_ids.join(
            df_term2gene,
            left_on="InputId",
            right_on="GeneId",
            how="inner",
        )
        .with_columns(pl.col("InputId").alias("GeneId"))
        .join(
            df_pathway.select("WikiPathwaysId", "PathwayName", "Species", "Url"),
            on="WikiPathwaysId",
            how="left",
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
