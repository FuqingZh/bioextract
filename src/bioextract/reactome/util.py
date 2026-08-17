from __future__ import annotations

import posixpath
import zipfile
from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from .constant import (
    COLS_MAPPING_RAW,
    COLS_PATHWAY_RAW,
    COLS_RELATION_RAW,
    SCHEMA_MAPPING_RAW,
    SCHEMA_PATHWAY_RAW,
    SCHEMA_RELATION_RAW,
    MappingRoleSpec,
)


def read_mapping_frame(file_uniprot2reactome: Path) -> pl.DataFrame:
    return read_mapping_family_frame(
        file_uniprot2reactome,
        columns=COLS_MAPPING_RAW,
        schema=SCHEMA_MAPPING_RAW,
        context="Reactome UniProt2Reactome file",
    )


def read_mapping_family_frame(
    file_mapping: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
    context: str,
) -> pl.DataFrame:
    """Read one official six-column Reactome mapping-family relation.

    The upstream files are literal TSV records, not CSV records. Quoting is
    therefore disabled so an embedded double quote remains part of the field,
    while exact-width and required-field checks keep malformed records out of
    canonical publication.
    """
    return _read_reactome_tsv(
        file_mapping,
        columns=columns,
        schema=schema,
        context=context,
        deduplicate=True,
    )


def scan_mapping_frame(file_uniprot2reactome: Path) -> pl.LazyFrame:
    return scan_mapping_family_frame(
        file_uniprot2reactome,
        columns=COLS_MAPPING_RAW,
        schema=SCHEMA_MAPPING_RAW,
        context="Reactome UniProt2Reactome file",
    )


def scan_mapping_family_frame(
    file_mapping: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
    context: str,
) -> pl.LazyFrame:
    """Return a strict-schema lazy scan for one mapping-family relation."""
    return _scan_reactome_tsv(
        file_mapping,
        columns=columns,
        schema=schema,
        context=context,
    )


def read_mapping_role_frame(
    file_mapping: Path,
    role: MappingRoleSpec,
) -> pl.DataFrame:
    """Read one registered Reactome mapping role."""
    return read_mapping_family_frame(
        file_mapping,
        columns=list(role.raw_columns),
        schema=dict.fromkeys(role.raw_columns, pl.String),
        context=f"Reactome {role.role} file",
    )


def scan_mapping_role_frame(
    file_mapping: Path,
    role: MappingRoleSpec,
) -> pl.LazyFrame:
    """Return a strict lazy scan for one registered mapping role."""
    return scan_mapping_family_frame(
        file_mapping,
        columns=list(role.raw_columns),
        schema=dict.fromkeys(role.raw_columns, pl.String),
        context=f"Reactome {role.role} file",
    )


def read_entity_pathway_frame(
    file_path: Path,
    *,
    source_columns: tuple[str, str, str],
    public_columns: tuple[str, str, str],
    context: str,
) -> pl.DataFrame:
    """Read one exact-header Reactome human entity relation."""
    lf = scan_entity_pathway_frame(
        file_path,
        source_columns=source_columns,
        public_columns=public_columns,
        context=context,
    )
    try:
        frame = lf.collect()
    except Exception as error:
        raise ValueError(f"{context} contains an invalid TSV record") from error
    _validate_required_values(frame, list(public_columns), context)
    return frame.unique(maintain_order=True).sort(list(public_columns))


def scan_entity_pathway_frame(
    file_path: Path,
    *,
    source_columns: tuple[str, str, str],
    public_columns: tuple[str, str, str],
    context: str,
) -> pl.LazyFrame:
    """Return a strict lazy scan for one headered human entity relation."""
    try:
        lf = pl.scan_csv(
            file_path,
            separator="\t",
            has_header=True,
            schema_overrides=dict.fromkeys(source_columns, pl.String),
            quote_char=None,
            truncate_ragged_lines=False,
        )
        observed_columns = lf.collect_schema().names()
    except Exception as error:
        raise ValueError(
            f"{context} must contain the exact ordered header: {list(source_columns)!r}"
        ) from error
    if observed_columns != list(source_columns):
        raise ValueError(
            f"{context} must contain the exact ordered header: {list(source_columns)!r}"
        )
    return lf.rename(dict(zip(source_columns, public_columns, strict=True)))


def read_gmt_frame(
    file_path: Path,
    *,
    public_columns: tuple[str, str, str],
    context: str,
) -> pl.DataFrame:
    """Read the one-member human Reactome GMT archive without extraction."""
    rows: list[tuple[str, str, str]] = []
    labels_by_pathway: dict[str, set[str]] = {}
    try:
        with zipfile.ZipFile(file_path) as archive:
            infos = archive.infolist()
            if len(infos) != 1:
                raise ValueError(f"{context} must contain exactly one file entry")
            info = infos[0]
            normalized_name = posixpath.normpath(info.filename)
            if (
                normalized_name != "ReactomePathways.gmt"
                or info.is_dir()
                or info.filename.startswith("/")
                or "\\" in info.filename
                or info.filename.startswith("../")
                or "/../" in info.filename
                or info.filename == ".."
                or info.flag_bits & 0x1
            ):
                raise ValueError(
                    f"{context} must contain one unencrypted ReactomePathways.gmt member"
                )
            with archive.open(info, "r") as binary_stream:
                for line_number, raw_line in enumerate(binary_stream, start=1):
                    try:
                        line = raw_line.decode("utf-8").rstrip("\r\n")
                    except UnicodeDecodeError as error:
                        raise ValueError(
                            f"{context} is not valid UTF-8 at line {line_number}"
                        ) from error
                    fields = line.split("\t")
                    if len(fields) < 3 or any(field == "" for field in fields):
                        raise ValueError(
                            f"{context} has an invalid GMT record at line {line_number}"
                        )
                    gene_set_name, pathway_id, *symbols = fields
                    labels_by_pathway.setdefault(pathway_id, set()).add(gene_set_name)
                    rows.extend(
                        (pathway_id, gene_set_name, symbol) for symbol in symbols
                    )
    except zipfile.BadZipFile as error:
        raise ValueError(f"{context} is not a valid ZIP archive") from error
    except (OSError, RuntimeError) as error:
        raise ValueError(f"{context} cannot be read") from error
    conflicting = sorted(
        pathway_id
        for pathway_id, labels in labels_by_pathway.items()
        if len(labels) > 1
    )
    if conflicting:
        raise ValueError(
            f"{context} has multiple gene-set labels for pathway {conflicting[0]}"
        )
    return (
        pl.DataFrame(
            {
                public_columns[0]: [row[0] for row in rows],
                public_columns[1]: [row[1] for row in rows],
                public_columns[2]: [row[2] for row in rows],
            },
            schema=dict.fromkeys(public_columns, pl.String),
        )
        .unique()
        .sort(list(public_columns))
    )


def read_pathway_frame(file_pathways: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_pathways,
        columns=COLS_PATHWAY_RAW,
        schema=SCHEMA_PATHWAY_RAW,
        context="Reactome pathways file",
    )


def scan_pathway_frame(file_pathways: Path) -> pl.LazyFrame:
    return _scan_reactome_tsv(
        file_pathways,
        columns=COLS_PATHWAY_RAW,
        schema=SCHEMA_PATHWAY_RAW,
        context="Reactome pathways file",
    )


def read_relation_frame(file_relations: Path) -> pl.DataFrame:
    return _read_reactome_tsv(
        file_relations,
        columns=COLS_RELATION_RAW,
        schema=SCHEMA_RELATION_RAW,
        context="Reactome pathway relations file",
    )


def scan_relation_frame(file_relations: Path) -> pl.LazyFrame:
    return _scan_reactome_tsv(
        file_relations,
        columns=COLS_RELATION_RAW,
        schema=SCHEMA_RELATION_RAW,
        context="Reactome pathway relations file",
    )


def filter_species_frame(df: pl.DataFrame, species: str | None) -> pl.DataFrame:
    if species is None:
        return df
    return df.filter(pl.col("species") == species)


def filter_relation_frame(
    df_relations: pl.DataFrame,
    df_pathways: pl.DataFrame,
) -> pl.DataFrame:
    df_pathway_ids = df_pathways.select("reactome_pathway_id").unique()
    return (
        df_relations.join(
            df_pathway_ids.rename(
                {"reactome_pathway_id": "parent_reactome_pathway_id"}
            ),
            on="parent_reactome_pathway_id",
            how="inner",
        )
        .join(
            df_pathway_ids.rename({"reactome_pathway_id": "child_reactome_pathway_id"}),
            on="child_reactome_pathway_id",
            how="inner",
        )
        .unique()
        .sort("parent_reactome_pathway_id", "child_reactome_pathway_id")
    )


def extract_mapping_frame(
    df_mapping: pl.DataFrame,
    df_input_ids: pl.DataFrame,
    *,
    df_group_membership: pl.DataFrame | None,
) -> pl.DataFrame:
    cols_out = [
        "input_id",
        "uniprot_id",
        "reactome_pathway_id",
        "pathway_name",
        "evidence_code",
        "species",
        "reactome_url",
    ]
    df_hits = (
        df_input_ids.join(
            df_mapping,
            left_on="input_id",
            right_on="uniprot_id",
            how="inner",
        )
        .with_columns(pl.col("input_id").alias("uniprot_id"))
        .select(cols_out)
        .unique()
        .sort(cols_out)
    )
    if df_group_membership is None:
        return df_hits
    grouped_cols = ["group_id", *cols_out]
    return (
        df_group_membership.join(df_hits, on="input_id", how="inner")
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


def extract_term2gene_frame(df_mapping: pl.DataFrame) -> pl.DataFrame:
    return (
        df_mapping.select("reactome_pathway_id", "uniprot_id")
        .unique()
        .sort("reactome_pathway_id", "uniprot_id")
    )


def extract_term2name_frame(df_pathways: pl.DataFrame) -> pl.DataFrame:
    return (
        df_pathways.select("reactome_pathway_id", "pathway_name", "species")
        .unique(subset=["reactome_pathway_id"])
        .sort("reactome_pathway_id")
    )


def _read_reactome_tsv(
    file_path: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
    context: str,
    deduplicate: bool = False,
) -> pl.DataFrame:
    lf = _scan_reactome_tsv(
        file_path,
        columns=columns,
        schema=schema,
        context=context,
    )
    try:
        df = lf.collect()
    except Exception as error:
        raise ValueError(f"{context} contains an invalid TSV record") from error
    _validate_required_values(df, columns, context)
    if deduplicate:
        df = df.unique(maintain_order=True)
    return df


def _scan_reactome_tsv(
    file_path: Path,
    *,
    columns: list[str],
    schema: SchemaDict,
    context: str,
) -> pl.LazyFrame:
    try:
        lf = pl.scan_csv(
            file_path,
            separator="\t",
            has_header=False,
            new_columns=columns,
            schema_overrides=schema,
            quote_char=None,
            truncate_ragged_lines=False,
        )
        observed_columns = lf.collect_schema().names()
    except Exception as error:
        raise ValueError(
            f"{context} must contain exactly {len(columns)} tab-separated fields"
        ) from error
    if observed_columns != columns:
        raise ValueError(
            f"{context} must contain exactly {len(columns)} tab-separated fields"
        )
    return lf


def _validate_required_values(
    df: pl.DataFrame,
    columns: list[str],
    context: str,
) -> None:
    if df.columns != columns:
        raise ValueError(
            f"{context} must contain exactly {len(columns)} tab-separated fields"
        )
    empty = df.select(
        pl.any_horizontal(
            [pl.col(column).is_null() | (pl.col(column) == "") for column in columns]
        ).any()
    ).item()
    if bool(empty):
        raise ValueError(f"{context} contains empty required fields")
