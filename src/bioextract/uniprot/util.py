from __future__ import annotations

import gzip
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO, TypedDict

import polars as pl

from bioextract._shared import RowWriter, create_tsv_writer, validate_required_cols

from .constant import (
    COLS_EGGNOG_XREF,
    COLS_IDMAPPING_SELECTED,
    COLS_SUBCELLULAR_LOCATION,
    REQUIRED_COLS_MAPPING,
    SCHEMA_EGGNOG_XREF,
    SCHEMA_MAPPING,
    SCHEMA_SUBCELLULAR_LOCATION,
)

_SUBCELLULAR_LOCATION_PREFIX = "CC   -!- SUBCELLULAR LOCATION:"
_CC_TOPIC_PREFIX = "CC   -!- "
_EVIDENCE_RE = re.compile(r"\{([^{}]+)\}")
_GENE_NAME_RE = re.compile(r"(?:^|;\s*)Name=([^;]+)")
_PROTEIN_FULL_RE = re.compile(r"RecName:\s+Full=([^;]+)")

type _EvidenceReference = tuple[str | None, str | None, str | None]


class _SubcellularLocationEntry(TypedDict):
    location: str | None
    note: str | None
    evidences: list[_EvidenceReference]


@dataclass(slots=True)
class _UniProtDatRecord:
    accessions: list[str] = field(default_factory=list)
    entry_name: str | None = None
    gene_name: str | None = None
    protein_name: str | None = None
    subcellular_location_comments: list[str] = field(default_factory=list)


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


def read_subcellular_location_frame(
    file_dat: Path,
    *,
    source_db: str,
) -> pl.DataFrame:
    """Extract a deduplicated accession-by-location-by-evidence table.

    Secondary accessions receive the same curated annotation as the primary
    accession. Locations without evidence retain one row with null evidence
    fields, keeping missing evidence distinct from a missing annotation.
    """
    rows: list[dict[str, str | None]] = []
    for record in iter_subcellular_location_records(file_dat):
        append_subcellular_location_rows(
            rows,
            record=record,
            source_db=source_db,
        )

    if not rows:
        return pl.DataFrame(schema=SCHEMA_SUBCELLULAR_LOCATION)
    return (
        pl.DataFrame(rows, schema=SCHEMA_SUBCELLULAR_LOCATION)
        .unique()
        .sort(COLS_SUBCELLULAR_LOCATION)
    )


def iter_subcellular_location_records(file_dat: Path) -> Iterable[_UniProtDatRecord]:
    """Yield UniProt flat-file records with identity and wrapped CC comments.

    Record termination uses ``//``. The final partial record is yielded only
    when it contains identity or subcellular-location content, which supports
    compact fixtures without weakening normal flat-file boundaries.
    """
    with open_uniprot_dat(file_dat) as handle:
        record = _UniProtDatRecord()
        is_subcellular_comment = False
        subcellular_comment_parts: list[str] = []
        for line in handle:
            if line.startswith("//"):
                if is_subcellular_comment:
                    record.subcellular_location_comments.append(
                        join_wrapped_comment_lines(subcellular_comment_parts)
                    )
                yield record
                record = _UniProtDatRecord()
                is_subcellular_comment = False
                subcellular_comment_parts = []
                continue

            if is_subcellular_comment and line.startswith(_CC_TOPIC_PREFIX):
                record.subcellular_location_comments.append(
                    join_wrapped_comment_lines(subcellular_comment_parts)
                )
                is_subcellular_comment = False
                subcellular_comment_parts = []

            if line.startswith(_SUBCELLULAR_LOCATION_PREFIX):
                is_subcellular_comment = True
                subcellular_comment_parts = [
                    line.removeprefix(_SUBCELLULAR_LOCATION_PREFIX).strip()
                ]
                continue
            if is_subcellular_comment and line.startswith("CC       "):
                subcellular_comment_parts.append(line[9:].strip())
                continue

            update_dat_record_identity(record, line)

        if is_subcellular_comment:
            record.subcellular_location_comments.append(
                join_wrapped_comment_lines(subcellular_comment_parts)
            )
        if (
            record.accessions
            or record.entry_name
            or record.subcellular_location_comments
        ):
            yield record


def update_dat_record_identity(record: _UniProtDatRecord, line: str) -> None:
    if line.startswith("ID   "):
        record.entry_name = line[5:].strip().split()[0]
    elif line.startswith("AC   "):
        record.accessions.extend(parse_accession_line(line))
    elif record.gene_name is None and line.startswith("GN   "):
        record.gene_name = parse_gene_name_line(line)
    elif record.protein_name is None and line.startswith("DE   RecName: Full="):
        record.protein_name = parse_recommended_protein_name_line(line)


def append_subcellular_location_rows(
    rows: list[dict[str, str | None]],
    *,
    record: _UniProtDatRecord,
    source_db: str,
) -> None:
    if not record.accessions or not record.subcellular_location_comments:
        return
    primary_accession = record.accessions[0]
    for comment in record.subcellular_location_comments:
        entries = parse_subcellular_location_comment(comment)
        for accession in record.accessions:
            for entry in entries:
                for evidence_code, evidence_source, evidence_id in entry["evidences"]:
                    rows.append(
                        {
                            "UniProtId": accession,
                            "PrimaryUniProtId": primary_accession,
                            "UniProtEntryName": record.entry_name,
                            "GeneName": record.gene_name,
                            "ProteinName": record.protein_name,
                            "SubcellularLocation": entry["location"],
                            "SubcellularLocationNote": entry["note"],
                            "EvidenceCode": evidence_code,
                            "EvidenceSource": evidence_source,
                            "EvidenceId": evidence_id,
                            "SourceDb": source_db,
                        }
                    )


def parse_subcellular_location_comment(
    comment: str,
) -> list[_SubcellularLocationEntry]:
    """Project one UniProt subcellular-location comment into long entries.

    Top-level location statements become separate entries; a trailing
    ``Note=`` applies to every entry. Each entry always contains at least one
    evidence tuple, using an all-null tuple when no evidence was supplied.
    """
    location_text, note_text = split_subcellular_location_note(comment)
    statements = split_top_level_periods(location_text)
    entries = [
        create_subcellular_location_entry(statement, note_text)
        for statement in statements
        if statement.strip()
    ]
    if not entries and note_text:
        entries.append(create_subcellular_location_entry("", note_text))
    return entries


def create_subcellular_location_entry(
    text: str,
    note_text: str | None,
) -> _SubcellularLocationEntry:
    location_text, evidences = extract_evidence_references(text)
    return {
        "location": location_text or None,
        "note": note_text,
        "evidences": evidences or [(None, None, None)],
    }


def split_subcellular_location_note(comment: str) -> tuple[str, str | None]:
    if "Note=" not in comment:
        return comment.strip(), None
    location_text, note_text = comment.split("Note=", 1)
    note_text, _ = extract_evidence_references(note_text)
    return location_text.strip(), note_text.strip().rstrip(".") or None


def join_wrapped_comment_lines(parts: list[str]) -> str:
    text = ""
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if not text:
            text = part
        elif text.endswith("-"):
            text += part
        else:
            text += f" {part}"
    return re.sub(r"\s+", " ", text).strip()


def split_top_level_periods(text: str) -> list[str]:
    """Split location statements on periods outside evidence braces."""
    statements: list[str] = []
    start = 0
    depth_brace = 0
    for index, character in enumerate(text):
        if character == "{":
            depth_brace += 1
        elif character == "}":
            depth_brace = max(0, depth_brace - 1)
        elif character == "." and depth_brace == 0:
            statement = text[start:index].strip()
            if statement:
                statements.append(statement)
            start = index + 1
    tail = text[start:].strip()
    if tail:
        statements.append(tail)
    return statements


def extract_evidence_references(
    text: str,
) -> tuple[str, list[_EvidenceReference]]:
    """Remove evidence blocks and return their parsed ECO/source references.

    Multiple comma-separated references in one brace block are preserved in
    source order. Unrecognized partial references retain the components that
    can be parsed rather than discarding the surrounding annotation text.
    """
    evidences: list[_EvidenceReference] = []
    for evidence_block in _EVIDENCE_RE.findall(text):
        for evidence_reference in evidence_block.split(","):
            evidence_reference = evidence_reference.strip()
            if not evidence_reference:
                continue
            evidence_code, evidence_source, evidence_id = parse_evidence_reference(
                evidence_reference
            )
            evidences.append((evidence_code, evidence_source, evidence_id))
    cleaned = _EVIDENCE_RE.sub("", text)
    cleaned = re.sub(r"\s+", " ", cleaned).strip(" ;")
    return cleaned, evidences


def parse_evidence_reference(
    evidence_reference: str,
) -> tuple[str | None, str | None, str | None]:
    evidence_code, separator, source_reference = evidence_reference.partition("|")
    evidence_code = evidence_code.strip() or None
    if not separator:
        return evidence_code, None, None
    evidence_source, source_separator, evidence_id = source_reference.strip().partition(
        ":"
    )
    if not source_separator:
        return evidence_code, evidence_source.strip() or None, None
    return evidence_code, evidence_source.strip() or None, evidence_id.strip() or None


def parse_gene_name_line(line: str) -> str | None:
    match = _GENE_NAME_RE.search(line[5:].strip())
    if match is None:
        return None
    gene_name, _ = extract_evidence_references(match.group(1))
    return gene_name or None


def parse_recommended_protein_name_line(line: str) -> str | None:
    match = _PROTEIN_FULL_RE.search(line[5:].strip())
    if match is None:
        return None
    protein_name, _ = extract_evidence_references(match.group(1))
    return protein_name or None


def write_eggnog_xref_tsv(
    file_dat: Path,
    path: Path,
    *,
    source_db: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle_out:
        writer = create_tsv_writer(handle_out)
        writer.writerow(COLS_EGGNOG_XREF)
        with open_uniprot_dat(file_dat) as handle_in:
            accessions: list[str] = []
            eggnog_xrefs: list[tuple[str, str]] = []
            for line in handle_in:
                if line.startswith("//"):
                    write_eggnog_xref_rows(
                        writer,
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
                writer,
                accessions=accessions,
                eggnog_xrefs=eggnog_xrefs,
                source_db=source_db,
            )


def write_subcellular_location_tsv(
    file_dat: Path,
    path: Path,
    *,
    source_db: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle_out:
        writer = create_tsv_writer(handle_out)
        writer.writerow(COLS_SUBCELLULAR_LOCATION)
        for record in iter_subcellular_location_records(file_dat):
            write_subcellular_location_rows(
                writer,
                record=record,
                source_db=source_db,
            )


def scan_eggnog_xref_tsv(file_tsv: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_tsv,
        separator="\t",
        has_header=True,
        schema_overrides=SCHEMA_EGGNOG_XREF,
    ).select(COLS_EGGNOG_XREF)


def scan_subcellular_location_tsv(file_tsv: Path) -> pl.LazyFrame:
    return (
        pl.scan_csv(
            file_tsv,
            separator="\t",
            has_header=True,
            schema_overrides=SCHEMA_SUBCELLULAR_LOCATION,
            infer_schema=False,
            null_values=[""],
        )
        .select(COLS_SUBCELLULAR_LOCATION)
        .unique()
        .sort(COLS_SUBCELLULAR_LOCATION)
    )


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
    writer: RowWriter,
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
            writer.writerow(
                [
                    accession,
                    primary_accession,
                    is_primary,
                    eggnog_og_id,
                    eggnog_level,
                    source_db,
                ]
            )


def write_subcellular_location_rows(
    writer: RowWriter,
    *,
    record: _UniProtDatRecord,
    source_db: str,
) -> None:
    rows: list[dict[str, str | None]] = []
    append_subcellular_location_rows(
        rows,
        record=record,
        source_db=source_db,
    )
    for row in rows:
        writer.writerow([row[column] or "" for column in COLS_SUBCELLULAR_LOCATION])


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
