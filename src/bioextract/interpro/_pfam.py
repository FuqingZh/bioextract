from __future__ import annotations

import gzip
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import polars as pl

from bioextract._tidy import (
    TidyAsset,
    TidyDataset,
    TidySource,
    calculate_file_sha256,
)

from .util import scan_protein2ipr_frame

PFAM_SCHEMA_VERSION = "interpro-pfam-v0.1"
PFAM_ID_PATTERN = re.compile(r"^PF[0-9]{5}$")

MEDIA_TYPE_TSV_GZIP = "application/gzip+tab-separated-values"
MEDIA_TYPE_XML_GZIP = "application/gzip+xml"

SCHEMA_PROTEIN_TERM = pl.Schema(
    {
        "UniProtId": pl.String,
        "PfamId": pl.String,
    }
)
SCHEMA_TERM = pl.Schema(
    {
        "PfamId": pl.String,
        "PfamName": pl.String,
    }
)
SCHEMA_TERM_XREF = pl.Schema(
    {
        "PfamId": pl.String,
        "InterProId": pl.String,
        "InterProName": pl.String,
        "InterProType": pl.String,
    }
)
_SCHEMA_XML_MAPPING = pl.Schema(
    {
        **SCHEMA_TERM_XREF,
        "PfamName": pl.String,
    }
)

_ASSETS = (
    TidyAsset(
        path="protein_term.parquet",
        kind="canonical",
        frame_name="protein_term",
    ),
    TidyAsset(path="term.parquet", kind="canonical", frame_name="term"),
    TidyAsset(
        path="term_xref.parquet",
        kind="canonical",
        frame_name="term_xref",
    ),
)


def build_pfam_tidy_dataset(
    *,
    file_protein2ipr: Path,
    file_interpro_xml: Path,
    include_source_hashes: bool = False,
) -> TidyDataset:
    """Build compact Pfam assets from explicitly assigned InterPro source roles.

    Args:
        file_protein2ipr: Protein-to-InterPro mapping source.
        file_interpro_xml: InterPro XML metadata source. Its official database
            metadata supplies the release identity.
        include_source_hashes: Whether to calculate source SHA-256 values for a
            later manifest write.

    Returns:
        A tidy dataset with lazy `protein_term`, `term`, and `term_xref` frames.
        Only Pfam IDs used by the protein mapping are retained.

    Raises:
        ValueError: If Pfam IDs or names are invalid, mapping relationships do
            not exist in the XML metadata, or a published frame would violate
            its schema.

    Notes:
        Paths declare logical roles only and never carry release identity.
        This boundary intentionally builds directly from the two raw files.
        Keep exact `InterProId + PfamId` validation here; deriving from a prior
        canonical parquet would add a publication prerequisite and could hide
        mismatched raw relationships. Pfam names come from XML member metadata,
        while InterPro entry names are retained only in `term_xref`.
    """
    interpro_version, df_xml_mapping = _read_pfam_xml_mapping(file_interpro_xml)
    _validate_xml_pfam_ids(df_xml_mapping)

    lf_raw_mapping = scan_protein2ipr_frame(file_protein2ipr)
    lf_raw_pfam = lf_raw_mapping.filter(
        pl.col("MemberDbId").str.starts_with("PF")
    ).select(
        pl.col("UniProtId").str.strip_chars().replace("", None),
        pl.col("MemberDbId").str.strip_chars().replace("", None).alias("PfamId"),
        pl.col("InterProId").str.strip_chars().replace("", None),
    )
    df_used_pair_index = (
        lf_raw_pfam.group_by("InterProId", "PfamId")
        .agg(
            pl.col("UniProtId").is_null().any().alias("HasMissingUniProtId"),
        )
        .collect(engine="streaming")
        .sort("InterProId", "PfamId", nulls_last=True)
    )
    _validate_raw_pfam_ids(df_used_pair_index)

    lf_used_pairs = df_used_pair_index.select("InterProId", "PfamId").lazy()
    lf_used_metadata = lf_used_pairs.join(
        df_xml_mapping.lazy(),
        on=["InterProId", "PfamId"],
        how="left",
    )
    _validate_used_xml_mapping(lf_used_metadata)

    lf_term = (
        lf_used_metadata.select(SCHEMA_TERM.names())
        .unique(keep="any", maintain_order=False)
        .sort("PfamId")
    )
    lf_term_xref = (
        lf_used_metadata.select(SCHEMA_TERM_XREF.names())
        .unique(keep="any", maintain_order=False)
        .sort("PfamId", "InterProId")
    )
    lf_protein_term = (
        lf_raw_pfam.join(
            lf_used_pairs,
            on=["InterProId", "PfamId"],
            how="semi",
        )
        .select(
            pl.col("UniProtId"),
            pl.col("PfamId"),
        )
        .unique(
            subset=SCHEMA_PROTEIN_TERM.names(),
            keep="any",
            maintain_order=False,
        )
    )
    frames = {
        "protein_term": lf_protein_term,
        "term": lf_term,
        "term_xref": lf_term_xref,
    }
    _validate_frame_schemas(frames)

    source_hashes = (
        {
            file_protein2ipr: calculate_file_sha256(file_protein2ipr),
            file_interpro_xml: calculate_file_sha256(file_interpro_xml),
        }
        if include_source_hashes
        else {}
    )
    return TidyDataset(
        resource_name="interpro",
        frames=frames,
        source=(
            TidySource(
                logical_name="protein_to_interpro",
                path=file_protein2ipr,
                media_type=MEDIA_TYPE_TSV_GZIP,
                sha256=source_hashes.get(file_protein2ipr),
            ),
            TidySource(
                logical_name="interpro_xml",
                path=file_interpro_xml,
                media_type=MEDIA_TYPE_XML_GZIP,
                sha256=source_hashes.get(file_interpro_xml),
            ),
        ),
        resource_schema_version=PFAM_SCHEMA_VERSION,
        source_schema_profile="interpro-pfam-bundle-v1",
        release_version=interpro_version,
        release_version_source="official_metadata",
        build_id_prefix=f"interpro-pfam-{interpro_version}",
        assets=_ASSETS,
    )


def _read_pfam_xml_mapping(file_interpro_xml: Path) -> tuple[str, pl.DataFrame]:
    interpro_versions: set[str] = set()
    rows: list[dict[str, str | None]] = []

    with gzip.open(file_interpro_xml, "rb") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag == "dbinfo":
                if elem.attrib.get("dbname") == "INTERPRO":
                    version = _clean_text(elem.attrib.get("version"))
                    if version is not None:
                        interpro_versions.add(version)
                elem.clear()
                continue
            if elem.tag != "interpro":
                continue

            interpro_id = _clean_text(elem.attrib.get("id"))
            interpro_type = _clean_text(elem.attrib.get("type"))
            elem_name = elem.find("name")
            interpro_name = (
                _clean_text("".join(elem_name.itertext()))
                if elem_name is not None
                else None
            )
            member_list = elem.find("member_list")
            if member_list is not None:
                for db_xref in member_list.findall("db_xref"):
                    if db_xref.attrib.get("db") != "PFAM":
                        continue
                    rows.append(
                        {
                            "PfamId": _clean_text(db_xref.attrib.get("dbkey")),
                            "PfamName": _clean_text(db_xref.attrib.get("name")),
                            "InterProId": interpro_id,
                            "InterProName": interpro_name,
                            "InterProType": interpro_type,
                        }
                    )
            elem.clear()

    if len(interpro_versions) != 1:
        raise ValueError(
            "InterPro XML must declare exactly one INTERPRO release version: "
            f"path={file_interpro_xml}, versions={sorted(interpro_versions)}"
        )
    if not rows:
        raise ValueError(
            f"InterPro XML contains no PFAM member references: {file_interpro_xml}"
        )

    return interpro_versions.pop(), pl.DataFrame(rows, schema=_SCHEMA_XML_MAPPING)


def read_interpro_release_version(file_interpro_xml: Path) -> str | None:
    versions: set[str] = set()
    with gzip.open(file_interpro_xml, "rb") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag == "dbinfo" and elem.attrib.get("dbname") == "INTERPRO":
                version = _clean_text(elem.attrib.get("version"))
                if version is not None:
                    versions.add(version)
            elem.clear()
    if not versions:
        return None
    if len(versions) != 1:
        raise ValueError(
            "InterPro XML must declare exactly one INTERPRO release version"
        )
    return versions.pop()


def _validate_xml_pfam_ids(df_xml_mapping: pl.DataFrame) -> None:
    invalid_ids = _invalid_pfam_ids(df_xml_mapping.get_column("PfamId").to_list())
    if invalid_ids:
        raise ValueError(
            f"InterPro XML contains invalid PFAM member IDs: {invalid_ids[:10]}"
        )


def _validate_raw_pfam_ids(df_pfam_id_stats: pl.DataFrame) -> None:
    if df_pfam_id_stats.is_empty():
        raise ValueError("InterPro protein2ipr contains no PFAM member IDs")

    invalid_ids = _invalid_pfam_ids(df_pfam_id_stats.get_column("PfamId").to_list())
    if invalid_ids:
        raise ValueError(
            f"InterPro protein2ipr contains invalid PFAM member IDs: {invalid_ids[:10]}"
        )
    missing_uniprot_ids = df_pfam_id_stats.filter("HasMissingUniProtId")
    if missing_uniprot_ids.height:
        raise ValueError(
            "InterPro PFAM rows contain missing UniProt IDs for: "
            f"{missing_uniprot_ids.select('InterProId', 'PfamId').head(10).rows()}"
        )


def _validate_used_xml_mapping(lf_used_metadata: pl.LazyFrame) -> None:
    cols_required = SCHEMA_TERM_XREF.names() + ["PfamName"]
    rows_incomplete = (
        lf_used_metadata.filter(pl.any_horizontal(pl.col(cols_required).is_null()))
        .select("InterProId", "PfamId")
        .unique(keep="any", maintain_order=False)
        .limit(10)
        .collect(engine="streaming")
    )
    if not rows_incomplete.is_empty():
        raise ValueError(
            "InterPro XML has incomplete PFAM metadata for mapped pairs: "
            f"{rows_incomplete.rows()}"
        )

    name_conflicts = (
        lf_used_metadata.group_by("PfamId")
        .agg(pl.col("PfamName").n_unique().alias("NameCount"))
        .filter(pl.col("NameCount") != 1)
        .sort("PfamId")
        .select("PfamId")
        .limit(10)
        .collect(engine="streaming")
    )
    if not name_conflicts.is_empty():
        raise ValueError(
            "InterPro XML assigns conflicting names to mapped PFAM IDs: "
            f"{name_conflicts.get_column('PfamId').to_list()}"
        )


def _validate_frame_schemas(frames: dict[str, pl.LazyFrame]) -> None:
    schema_by_frame = {
        "protein_term": SCHEMA_PROTEIN_TERM,
        "term": SCHEMA_TERM,
        "term_xref": SCHEMA_TERM_XREF,
    }
    for frame_name, schema_expected in schema_by_frame.items():
        schema_actual = frames[frame_name].collect_schema()
        if schema_actual != schema_expected:
            raise ValueError(
                "InterPro PFAM frame schema mismatch: "
                f"frame={frame_name!r}, expected={schema_expected!r}, "
                f"actual={schema_actual!r}"
            )


def _invalid_pfam_ids(values: list[str | None]) -> list[str | None]:
    return [
        value
        for value in values
        if not isinstance(value, str) or PFAM_ID_PATTERN.fullmatch(value) is None
    ]


def _clean_text(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    return value or None
