from __future__ import annotations

import bz2
import csv
import gzip
import hashlib
import lzma
import re
import sqlite3
import tempfile
import zlib
from collections.abc import Generator, Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

import polars as pl
from polars._typing import SchemaDict

from bioextract._publication import (
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    preflight_publication_destination,
    write_duckdb_publication,
)

from .util import join_wrapped_comment_lines, parse_subcellular_location_comment

RESOURCE_SCHEMA_VERSION = "uniprot-knowledgebase-duckdb-v1"
SOURCE_SCHEMA_PROFILE = "uniprotkb-flat-file-v1"

_EVIDENCE = re.compile(r"\{([^{}]+)\}")
_ISO_ID = re.compile(r"IsoId=([^;]+);")
_ISO_SEQUENCE = re.compile(r"Sequence=([^;]+);")
_VSP_IDENTIFIER = re.compile(r"VSP_[0-9]{6,10}")
_VSP_ID = re.compile(r"\b(VSP_[0-9]{6,10})\b")
_VSP_SEQUENCE = re.compile(r"VSP_[0-9]{6,10}(?:,\s*VSP_[0-9]{6,10})*")
_ISOFORM_SCOPE = re.compile(r"^(.*?)\s+\[([A-Z0-9-]+)\]$")
_CRC64_POLYNOMIAL = 0xD800000000000000
_MOLECULAR_WEIGHT_WATER = 180153
_RESIDUE_MOLECULAR_WEIGHT = {
    "A": 710788,
    "B": 1146532,
    "C": 1031388,
    "D": 1150886,
    "E": 1291155,
    "F": 1471766,
    "G": 570519,
    "H": 1371411,
    "I": 1131594,
    "K": 1281741,
    "L": 1131594,
    "M": 1311926,
    "N": 1141038,
    "O": 2373000,
    "P": 971167,
    "Q": 1281307,
    "R": 1561875,
    "S": 870782,
    "T": 1011051,
    "U": 1500400,
    "V": 991326,
    "W": 1862132,
    "X": 1113306,
    "Y": 1631760,
    "Z": 1287473,
}
_LEGACY_EXPASY_MOLECULAR_WEIGHT = {
    "O": 2373018,
    "U": 1500388,
}


def _crc64_table() -> tuple[int, ...]:
    values: list[int] = []
    for byte in range(256):
        value = byte
        for _ in range(8):
            value = (value >> 1) ^ _CRC64_POLYNOMIAL if value & 1 else value >> 1
        values.append(value)
    return tuple(values)


_CRC64_TABLE = _crc64_table()


@dataclass(slots=True)
class _Record:
    number: int
    entry_name: str | None = None
    reviewed: bool | None = None
    sequence_length: int | None = None
    molecular_weight: int | None = None
    crc64: str | None = None
    accessions: list[str] = field(default_factory=list[str])
    taxon_id: str | None = None
    protein_existence: str | None = None
    sequence_version: int | None = None
    entry_version: int | None = None
    sequence: list[str] = field(default_factory=list[str])
    de_parts: list[str] = field(default_factory=list[str])
    gn_parts: list[str] = field(default_factory=list[str])
    names: list[tuple[str, str]] = field(default_factory=list[tuple[str, str]])
    gene_names: list[tuple[str, str]] = field(default_factory=list[tuple[str, str]])
    ec_numbers: list[str] = field(default_factory=list[str])
    cross_references: list[tuple[str, list[str]]] = field(
        default_factory=list[tuple[str, list[str]]]
    )
    comments: list[tuple[str, str]] = field(default_factory=list[tuple[str, str]])
    keywords: list[str] = field(default_factory=list[str])
    isoforms: list[tuple[list[str], str | None, str]] = field(
        default_factory=list[tuple[list[str], str | None, str]]
    )
    variations: list[tuple[str, int | None, int | None, str]] = field(
        default_factory=list[tuple[str, int | None, int | None, str]]
    )
    molecular_weight_profile: str = "shared"


TABLE_SCHEMAS: Mapping[str, SchemaDict] = {
    "protein": {
        "primary_accession": pl.String,
        "entry_name": pl.String,
        "is_reviewed": pl.Boolean,
        "taxon_id": pl.String,
        "protein_existence": pl.String,
        "sequence_length": pl.Int64,
        "molecular_weight": pl.Int64,
        "sequence_version": pl.Int64,
        "entry_version": pl.Int64,
    },
    "protein_accession": {
        "primary_accession": pl.String,
        "accession": pl.String,
        "accession_order": pl.Int64,
        "is_primary": pl.Boolean,
    },
    "protein_sequence": {
        "sequence_id": pl.String,
        "primary_accession": pl.String,
        "sequence_type": pl.String,
        "sequence": pl.String,
        "length": pl.Int64,
        "crc64": pl.String,
        "sha256": pl.String,
    },
    "protein_isoform": {
        "isoform_id": pl.String,
        "primary_accession": pl.String,
        "name": pl.String,
        "isoform_order": pl.Int64,
        "sequence_status": pl.String,
        "sequence_id": pl.String,
    },
    "protein_isoform_identifier": {
        "primary_accession": pl.String,
        "isoform_id": pl.String,
        "identifier": pl.String,
        "identifier_order": pl.Int64,
        "is_main": pl.Boolean,
    },
    "protein_sequence_variation": {
        "variation_id": pl.String,
        "primary_accession": pl.String,
        "start_position": pl.Int64,
        "end_position": pl.Int64,
        "note": pl.String,
    },
    "protein_isoform_variation": {
        "primary_accession": pl.String,
        "isoform_id": pl.String,
        "variation_id": pl.String,
        "variation_order": pl.Int64,
    },
    "protein_name": {
        "primary_accession": pl.String,
        "name_type": pl.String,
        "name": pl.String,
        "name_order": pl.Int64,
    },
    "gene_name": {
        "primary_accession": pl.String,
        "name_type": pl.String,
        "name": pl.String,
        "name_order": pl.Int64,
    },
    "protein_ec_number": {
        "primary_accession": pl.String,
        "ec_number": pl.String,
    },
    "protein_go_annotation": {
        "primary_accession": pl.String,
        "go_id": pl.String,
        "aspect": pl.String,
        "term_name": pl.String,
        "evidence_code": pl.String,
        "evidence_source": pl.String,
    },
    "protein_cross_reference": {
        "primary_accession": pl.String,
        "database": pl.String,
        "external_id": pl.String,
        "properties": pl.String,
        "isoform_id": pl.String,
    },
    "protein_comment": {
        "comment_id": pl.String,
        "primary_accession": pl.String,
        "comment_type": pl.String,
        "comment_text": pl.String,
    },
    "protein_subcellular_location": {
        "comment_id": pl.String,
        "primary_accession": pl.String,
        "location": pl.String,
        "note": pl.String,
    },
    "protein_keyword": {
        "primary_accession": pl.String,
        "keyword": pl.String,
        "keyword_order": pl.Int64,
    },
    "protein_identifier": {
        "primary_accession": pl.String,
        "namespace": pl.String,
        "identifier": pl.String,
    },
}


def write_knowledgebase(
    *,
    entries: Path,
    canonical_sequences: Path | None,
    isoform_sequences: Path | None,
    release_version: str | None,
    path: Path,
    if_exists: str,
) -> DuckDBWriteResult:
    preflight_publication_destination(path, if_exists=if_exists)
    with tempfile.TemporaryDirectory(prefix="bioextract-uniprot-kb-") as temp_dir:
        spool_dir = Path(temp_dir)
        writers, handles = _open_spools(spool_dir)
        database_index = sqlite3.connect(spool_dir / "validation.sqlite")
        _create_validation_index(database_index)
        molecular_weight_profile: str | None = None
        try:
            for record in _iter_records(entries):
                if record.molecular_weight_profile != "shared":
                    if (
                        molecular_weight_profile is not None
                        and molecular_weight_profile != record.molecular_weight_profile
                    ):
                        raise ValueError(
                            "UniProtKB records require conflicting molecular-weight "
                            "models: current-core and legacy-expasy"
                        )
                    molecular_weight_profile = record.molecular_weight_profile
                _write_record(
                    record,
                    writers,
                    database_index=database_index,
                )
            database_index.commit()
            if canonical_sequences is not None:
                _validate_canonical_fasta(canonical_sequences, database_index)
            if isoform_sequences is not None:
                _write_isoform_fasta(
                    isoform_sequences,
                    writers["protein_sequence"],
                    database_index,
                )
            _write_isoform_relations(writers, database_index)
        finally:
            database_index.close()
            for handle in handles:
                handle.close()

        relations = tuple(
            RelationSpec(
                table_name=table,
                frame=pl.scan_csv(
                    spool_dir / f"{table}.tsv",
                    separator="\t",
                    schema_overrides=schema,
                    infer_schema=False,
                    null_values=[""],
                ),
            )
            for table, schema in TABLE_SCHEMAS.items()
        )
        sources = [
            SourceFileRecord("knowledgebase_entries", entries, _media_type(entries))
        ]
        if canonical_sequences is not None:
            sources.append(
                SourceFileRecord(
                    "canonical_fasta",
                    canonical_sequences,
                    _media_type(canonical_sequences),
                )
            )
        if isoform_sequences is not None:
            sources.append(
                SourceFileRecord(
                    "isoform_fasta",
                    isoform_sequences,
                    _media_type(isoform_sequences),
                )
            )
        return write_duckdb_publication(
            relations,
            path,
            resource_name="uniprot",
            resource_schema_version=RESOURCE_SCHEMA_VERSION,
            source_schema_profile=SOURCE_SCHEMA_PROFILE,
            sources=sources,
            scope="swiss-prot",
            release_version=release_version,
            release_version_source="caller" if release_version is not None else None,
            extra_metadata={
                "bioextract.source_schema_validation": "passed",
                "bioextract.molecular_weight_validation_model": (
                    molecular_weight_profile or "compatible-current-and-legacy"
                ),
                "bioextract.source_profile.knowledgebase_entries": (
                    "uniprotkb-swiss-prot-dat-v1"
                ),
                "bioextract.source_profile.canonical_fasta": (
                    "uniprotkb-canonical-fasta-v1"
                    if canonical_sequences is not None
                    else "absent"
                ),
                "bioextract.source_profile.isoform_fasta": (
                    "uniprotkb-varsplic-fasta-v1"
                    if isoform_sequences is not None
                    else "absent"
                ),
                "bioextract.capability.canonical_sequences": "true",
                "bioextract.capability.isoform_definitions": "true",
                "bioextract.capability.isoform_sequences": str(
                    isoform_sequences is not None
                ).lower(),
                "bioextract.canonical_fasta_validated": str(
                    canonical_sequences is not None
                ).lower(),
            },
            if_exists=if_exists,
        )


def _create_validation_index(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
            CREATE TABLE canonical(accession TEXT PRIMARY KEY, sequence TEXT NOT NULL);
            CREATE TABLE accession(
                accession TEXT NOT NULL,
                primary_accession TEXT NOT NULL,
                PRIMARY KEY (primary_accession, accession)
            );
            CREATE TABLE isoform(
                isoform_id TEXT NOT NULL,
                primary_accession TEXT NOT NULL,
                name TEXT,
                isoform_order INTEGER NOT NULL,
                sequence_status TEXT NOT NULL,
                sequence_id TEXT,
                PRIMARY KEY (primary_accession, isoform_id)
            );
            CREATE TABLE isoform_identifier(
                primary_accession TEXT NOT NULL,
                isoform_id TEXT NOT NULL,
                identifier TEXT NOT NULL,
                identifier_order INTEGER NOT NULL,
                is_main INTEGER NOT NULL,
                PRIMARY KEY (primary_accession, isoform_id, identifier),
                UNIQUE (primary_accession, identifier)
            );
            CREATE INDEX isoform_identifier_lookup
                ON isoform_identifier(identifier);
            CREATE TABLE isoform_variation(
                primary_accession TEXT NOT NULL,
                isoform_id TEXT NOT NULL,
                variation_id TEXT NOT NULL,
                variation_order INTEGER NOT NULL,
                PRIMARY KEY (primary_accession, isoform_id, variation_id)
            );
            CREATE TABLE variation(variation_id TEXT PRIMARY KEY);
            CREATE TABLE seen_isoform_sequence(
                primary_accession TEXT NOT NULL,
                isoform_id TEXT NOT NULL,
                PRIMARY KEY (primary_accession, isoform_id)
            );
        """
    )


def _open_spools(
    directory: Path,
) -> tuple[dict[str, csv.DictWriter[str]], list[IO[str]]]:
    writers: dict[str, csv.DictWriter[str]] = {}
    handles: list[IO[str]] = []
    for table, schema in TABLE_SCHEMAS.items():
        handle = (directory / f"{table}.tsv").open("w", encoding="utf-8", newline="")
        writer = csv.DictWriter(handle, fieldnames=list(schema), delimiter="\t")
        writer.writeheader()
        writers[table] = writer
        handles.append(handle)
    return writers, handles


def _iter_records(path: Path) -> Iterator[_Record]:
    with _open_text(path) as handle:
        record = _Record(number=1)
        in_sequence = False
        comment_type: str | None = None
        comment_parts: list[str] = []
        isoform_block = ""
        variation: tuple[str, int | None, int | None, str] | None = None
        variation_qualifier: str | None = None
        saw_record = False
        for line in handle:
            if line.startswith("//"):
                _flush_comment(record, comment_type, comment_parts)
                _parse_isoforms(record, isoform_block)
                _parse_de(record, " ".join(record.de_parts))
                _parse_gn(record, " ".join(record.gn_parts))
                if variation is not None:
                    record.variations.append(variation)
                _validate_record(record)
                yield record
                saw_record = True
                record = _Record(number=record.number + 1)
                in_sequence = False
                comment_type = None
                comment_parts = []
                isoform_block = ""
                variation = None
                variation_qualifier = None
                continue
            if in_sequence and line.startswith("     "):
                record.sequence.append("".join(line.split()))
                continue
            in_sequence = False
            if comment_type is not None and not line.startswith("CC   "):
                _flush_comment(record, comment_type, comment_parts)
                comment_type = None
                comment_parts = []
            feature_key = line[5:13].strip() if line.startswith("FT   ") else ""
            if feature_key and feature_key != "VAR_SEQ" and variation is not None:
                record.variations.append(variation)
                variation = None
                variation_qualifier = None
            if line.startswith("ID   "):
                parts = line[5:].split()
                if len(parts) < 3:
                    raise ValueError(
                        f"Invalid UniProt ID line in record {record.number}"
                    )
                record.entry_name = parts[0]
                record.reviewed = parts[1].rstrip(";") == "Reviewed"
                record.sequence_length = int(parts[-2])
            elif line.startswith("AC   "):
                record.accessions.extend(
                    value.strip() for value in line[5:].split(";") if value.strip()
                )
            elif line.startswith("OX   "):
                match = re.search(r"NCBI_TaxID=([0-9]+)", line)
                if match:
                    record.taxon_id = match.group(1)
            elif line.startswith("PE   "):
                record.protein_existence = line[5:].strip().rstrip(".;")
            elif line.startswith("DT   ") and "sequence version" in line:
                match = re.search(r"sequence version ([0-9]+)", line)
                record.sequence_version = int(match.group(1)) if match else None
            elif line.startswith("DT   ") and "entry version" in line:
                match = re.search(r"entry version ([0-9]+)", line)
                record.entry_version = int(match.group(1)) if match else None
            elif line.startswith("DE   "):
                record.de_parts.append(line[5:].strip())
            elif line.startswith("GN   "):
                record.gn_parts.append(line[5:].strip())
            elif line.startswith("DR   "):
                parts = [part.strip().rstrip(".") for part in line[5:].split(";")]
                if len(parts) >= 2:
                    record.cross_references.append((parts[0], parts[1:]))
            elif line.startswith("KW   "):
                record.keywords.extend(
                    value.strip().rstrip(".")
                    for value in line[5:].split(";")
                    if value.strip()
                )
            elif line.startswith("CC   -!- "):
                _flush_comment(record, comment_type, comment_parts)
                payload = line.removeprefix("CC   -!- ")
                comment_type, _, first = payload.partition(":")
                comment_parts = [first.strip()]
            elif line.startswith("CC       ") and comment_type is not None:
                comment_parts.append(line[9:].strip())
            elif line.startswith("SQ"):
                match = re.fullmatch(
                    r"SQ {3}SEQUENCE +([0-9]+) AA; +([0-9]+) MW; +"
                    r"([0-9A-F]{16}) CRC64;",
                    line.rstrip("\r\n"),
                )
                if match is None:
                    raise ValueError(
                        f"Invalid UniProt SQ line in record {record.number}"
                    )
                if int(match.group(1)) != record.sequence_length:
                    raise ValueError(
                        f"ID/SQ sequence length mismatch in record {record.number}"
                    )
                record.molecular_weight = int(match.group(2))
                record.crc64 = match.group(3)
                in_sequence = True
            elif line.startswith("FT   VAR_SEQ"):
                if variation is not None:
                    record.variations.append(variation)
                positions = re.findall(r"[0-9]+", line[13:35])
                variation = (
                    "",
                    int(positions[0]) if positions else None,
                    int(positions[-1]) if positions else None,
                    "",
                )
                variation_qualifier = None
            elif line.startswith("FT") and "/note=" in line and variation is not None:
                raw_note = line.split("=", 1)[1].strip()
                note = raw_note.removeprefix('"').removesuffix('"')
                variation = (
                    variation[0],
                    variation[1],
                    variation[2],
                    note,
                )
                variation_qualifier = None if raw_note.endswith('"') else "note"
            elif line.startswith("FT") and "/id=" in line and variation is not None:
                variation_id = line.split("=", 1)[1].strip().strip('"')
                if variation[0]:
                    raise ValueError(
                        f"Duplicate UniProt VAR_SEQ /id in record {record.number}"
                    )
                if _VSP_IDENTIFIER.fullmatch(variation_id) is None:
                    raise ValueError(
                        "Invalid UniProt VAR_SEQ /id in record "
                        f"{record.number}: {variation_id}"
                    )
                variation = (
                    variation_id,
                    variation[1],
                    variation[2],
                    variation[3],
                )
                variation_qualifier = None
            elif (
                line.startswith("FT") and "/evidence=" in line and variation is not None
            ):
                variation_qualifier = "evidence"
            elif (
                line.startswith("FT                   ")
                and variation is not None
                and variation[3]
                and variation_qualifier == "note"
            ):
                continuation = line[21:].strip()
                variation = (
                    variation[0],
                    variation[1],
                    variation[2],
                    f"{variation[3]} {continuation.removesuffix('"')}",
                )
                if continuation.endswith('"'):
                    variation_qualifier = None
            if comment_type == "ALTERNATIVE PRODUCTS" and line.startswith("CC   "):
                isoform_block += f" {line[9:].strip()}"
        if record.entry_name is not None:
            raise ValueError("Unterminated UniProtKB flat-file record")
        if not saw_record:
            raise ValueError("Declared entries input contains no UniProtKB records")


def _validate_record(record: _Record) -> None:
    missing: list[str] = []
    if record.entry_name is None:
        missing.append("ID")
    if record.reviewed is False:
        raise ValueError(
            f"Swiss-Prot source profile rejects Unreviewed record {record.number}"
        )
    if not record.accessions:
        missing.append("AC")
    if record.taxon_id is None:
        missing.append("OX")
    sequence = "".join(record.sequence)
    if record.sequence_length is None or not sequence:
        missing.append("SQ")
    if missing:
        raise ValueError(
            f"UniProtKB record {record.number} is missing required fields: {missing}"
        )
    if len(set(record.accessions)) != len(record.accessions):
        raise ValueError(f"Duplicate accession in UniProtKB record {record.number}")
    if any(residue not in _RESIDUE_MOLECULAR_WEIGHT for residue in sequence):
        raise ValueError(
            f"Invalid UniProt sequence characters in record {record.number}"
        )
    if len(sequence) != record.sequence_length:
        raise ValueError(
            f"Sequence length mismatch in UniProtKB record {record.number}"
        )
    actual_crc64 = _calculate_crc64(sequence)
    if record.crc64 != actual_crc64:
        raise ValueError(
            f"UniProt CRC64 mismatch in record {record.number}: "
            f"expected={record.crc64}, actual={actual_crc64}"
        )
    current_weight = _calculate_molecular_weight(sequence)
    legacy_weight = _calculate_molecular_weight(sequence, legacy_expasy=True)
    matches_current = record.molecular_weight == current_weight
    matches_legacy = record.molecular_weight == legacy_weight
    if not matches_current and not matches_legacy:
        raise ValueError(
            f"UniProt molecular weight mismatch in record {record.number}: "
            f"expected={record.molecular_weight}, "
            f"current={current_weight}, legacy_expasy={legacy_weight}"
        )
    if matches_current != matches_legacy:
        record.molecular_weight_profile = (
            "current-core" if matches_current else "legacy-expasy"
        )
    if any(not variation_id for variation_id, _start, _end, _note in record.variations):
        raise ValueError(
            f"UniProt VAR_SEQ feature is missing /id in record {record.number}"
        )


def _write_record(
    record: _Record,
    writers: Mapping[str, csv.DictWriter[str]],
    *,
    database_index: sqlite3.Connection,
) -> None:
    primary = record.accessions[0]
    sequence = "".join(record.sequence)
    try:
        database_index.execute(
            "INSERT INTO canonical VALUES (?, ?)", (primary, sequence)
        )
    except sqlite3.IntegrityError as error:
        raise ValueError(f"Duplicate primary UniProt accession: {primary}") from error
    writers["protein"].writerow(
        {
            "primary_accession": primary,
            "entry_name": record.entry_name,
            "is_reviewed": record.reviewed,
            "taxon_id": record.taxon_id,
            "protein_existence": record.protein_existence,
            "sequence_length": len(sequence),
            "molecular_weight": record.molecular_weight,
            "sequence_version": record.sequence_version,
            "entry_version": record.entry_version,
        }
    )
    sequence_id = f"{primary}:canonical"
    writers["protein_sequence"].writerow(
        {
            "sequence_id": sequence_id,
            "primary_accession": primary,
            "sequence_type": "canonical",
            "sequence": sequence,
            "length": len(sequence),
            "crc64": record.crc64,
            "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
    )
    for index, accession in enumerate(record.accessions):
        try:
            database_index.execute(
                "INSERT INTO accession VALUES (?, ?)", (accession, primary)
            )
        except sqlite3.IntegrityError as error:
            raise ValueError(
                f"Duplicate UniProt accession within record {record.number}: "
                f"{accession}"
            ) from error
        writers["protein_accession"].writerow(
            {
                "primary_accession": primary,
                "accession": accession,
                "accession_order": index + 1,
                "is_primary": index == 0,
            }
        )
        writers["protein_identifier"].writerow(
            {
                "primary_accession": primary,
                "namespace": "uniprot",
                "identifier": accession,
            }
        )
    writers["protein_identifier"].writerow(
        {
            "primary_accession": primary,
            "namespace": "entry_name",
            "identifier": record.entry_name,
        }
    )
    for index, (name_type, name) in enumerate(record.names, start=1):
        writers["protein_name"].writerow(
            {
                "primary_accession": primary,
                "name_type": name_type,
                "name": name,
                "name_order": index,
            }
        )
    for ec_number in record.ec_numbers:
        writers["protein_ec_number"].writerow(
            {"primary_accession": primary, "ec_number": ec_number}
        )
    for index, (name_type, name) in enumerate(record.gene_names, start=1):
        writers["gene_name"].writerow(
            {
                "primary_accession": primary,
                "name_type": name_type,
                "name": name,
                "name_order": index,
            }
        )
        writers["protein_identifier"].writerow(
            {
                "primary_accession": primary,
                "namespace": "gene_name",
                "identifier": name,
            }
        )
    official_isoform_ids = {
        identifier
        for identifiers, _name, _status in record.isoforms
        for identifier in identifiers
    }
    for database, fields in record.cross_references:
        if database == "GO":
            writers["protein_go_annotation"].writerow(
                {
                    "primary_accession": primary,
                    "go_id": fields[0],
                    "aspect": fields[1].split(":", 1)[0] if len(fields) > 1 else None,
                    "term_name": fields[1].split(":", 1)[-1]
                    if len(fields) > 1
                    else None,
                    "evidence_code": fields[2].split(":", 1)[0]
                    if len(fields) > 2
                    else None,
                    "evidence_source": fields[2].split(":", 1)[-1]
                    if len(fields) > 2 and ":" in fields[2]
                    else None,
                }
            )
            continue
        relation_fields = list(fields)
        isoform_scope = None
        if relation_fields:
            scope_match = _ISOFORM_SCOPE.fullmatch(relation_fields[-1])
            if scope_match is not None and scope_match.group(2) in official_isoform_ids:
                relation_fields[-1] = scope_match.group(1).strip().rstrip(".")
                isoform_scope = scope_match.group(2)
        writers["protein_cross_reference"].writerow(
            {
                "primary_accession": primary,
                "database": database,
                "external_id": relation_fields[0],
                "properties": "; ".join(relation_fields[1:]) or None,
                "isoform_id": isoform_scope,
            }
        )
        namespace = {
            "GeneID": "gene_id",
            "RefSeq": "refseq",
            "Ensembl": "ensembl",
        }.get(database)
        if namespace:
            for identifier in relation_fields:
                if identifier and identifier != "-":
                    writers["protein_identifier"].writerow(
                        {
                            "primary_accession": primary,
                            "namespace": namespace,
                            "identifier": identifier,
                        }
                    )
    for index, (comment_type, text) in enumerate(record.comments, start=1):
        comment_id = f"{primary}:comment:{index}"
        writers["protein_comment"].writerow(
            {
                "comment_id": comment_id,
                "primary_accession": primary,
                "comment_type": comment_type,
                "comment_text": text,
            }
        )
        if comment_type == "SUBCELLULAR LOCATION":
            for location in parse_subcellular_location_comment(text):
                writers["protein_subcellular_location"].writerow(
                    {
                        "comment_id": comment_id,
                        "primary_accession": primary,
                        "location": location["location"],
                        "note": location["note"],
                    }
                )
    for index, keyword in enumerate(record.keywords, start=1):
        writers["protein_keyword"].writerow(
            {
                "primary_accession": primary,
                "keyword": keyword,
                "keyword_order": index,
            }
        )
    for isoform_order, (identifiers, name, status) in enumerate(
        record.isoforms, start=1
    ):
        isoform_id = identifiers[0]
        variation_ids = _VSP_ID.findall(status)
        normalized_status = (
            status
            if status in {"Displayed", "External", "Not described"}
            else "Alternative"
        )
        try:
            database_index.execute(
                "INSERT INTO isoform VALUES (?, ?, ?, ?, ?, ?)",
                (
                    isoform_id,
                    primary,
                    name,
                    isoform_order,
                    normalized_status,
                    sequence_id if normalized_status == "Displayed" else None,
                ),
            )
        except sqlite3.IntegrityError as error:
            raise ValueError(
                f"Duplicate UniProt isoform product in one entry: {primary}/{isoform_id}"
            ) from error
        for identifier_order, identifier in enumerate(identifiers, start=1):
            try:
                database_index.execute(
                    "INSERT INTO isoform_identifier VALUES (?, ?, ?, ?, ?)",
                    (
                        primary,
                        isoform_id,
                        identifier,
                        identifier_order,
                        identifier_order == 1,
                    ),
                )
            except sqlite3.IntegrityError as error:
                raise ValueError(
                    "Duplicate UniProt isoform identifier in one entry: "
                    f"{primary}/{identifier}"
                ) from error
            writers["protein_identifier"].writerow(
                {
                    "primary_accession": primary,
                    "namespace": "isoform_id",
                    "identifier": identifier,
                }
            )
        for variation_order, variation_id in enumerate(variation_ids, start=1):
            database_index.execute(
                "INSERT INTO isoform_variation VALUES (?, ?, ?, ?)",
                (primary, isoform_id, variation_id, variation_order),
            )
    for variation_id, start, end, note in record.variations:
        try:
            database_index.execute("INSERT INTO variation VALUES (?)", (variation_id,))
        except sqlite3.IntegrityError as error:
            raise ValueError(
                f"Duplicate UniProt VAR_SEQ identifier: {variation_id}"
            ) from error
        writers["protein_sequence_variation"].writerow(
            {
                "variation_id": variation_id,
                "primary_accession": primary,
                "start_position": start,
                "end_position": end,
                "note": note,
            }
        )


def _parse_de(record: _Record, payload: str) -> None:
    record.ec_numbers.extend(re.findall(r"(?:^|; )EC=([0-9n.-]+)", payload))
    name_type_by_source = {
        "RecName": "recommended",
        "AltName": "alternative",
        "SubName": "submitted",
    }
    for match in re.finditer(r"(RecName|AltName|SubName):\s+Full=([^;]+)", payload):
        record.names.append(
            (
                name_type_by_source[match.group(1)],
                _EVIDENCE.sub("", match.group(2)).strip(),
            )
        )


def _parse_gn(record: _Record, payload: str) -> None:
    normalized = payload.replace("; and ", "; ")
    name_type_by_key = {
        "Name": "primary",
        "Synonyms": "synonym",
        "OrderedLocusNames": "ordered_locus",
        "ORFNames": "orf",
    }
    for match in re.finditer(
        r"(?:^|;\s+)(Name|Synonyms|OrderedLocusNames|ORFNames)=([^;]+)",
        normalized,
    ):
        clean_values = _EVIDENCE.sub("", match.group(2))
        record.gene_names.extend(
            (
                name_type_by_key[match.group(1)],
                value.strip(),
            )
            for value in clean_values.split(",")
            if value.strip()
        )


def _flush_comment(
    record: _Record, comment_type: str | None, comment_parts: list[str]
) -> None:
    if comment_type is not None:
        record.comments.append(
            (comment_type, join_wrapped_comment_lines(comment_parts))
        )


def _parse_isoforms(record: _Record, text: str) -> None:
    for block in text.split("Name=")[1:]:
        name, _, rest = block.partition(";")
        iso_ids = _ISO_ID.search(rest)
        sequence = _ISO_SEQUENCE.search(rest)
        if iso_ids is None or sequence is None:
            missing = "IsoId" if iso_ids is None else "Sequence"
            raise ValueError(
                "UniProt alternative-products Name block is missing "
                f"{missing} in record {record.number}"
            )
        if iso_ids and sequence:
            sequence_value = sequence.group(1).strip()
            if sequence_value not in {"Displayed", "External", "Not described"} and (
                _VSP_SEQUENCE.fullmatch(sequence_value) is None
            ):
                raise ValueError(
                    "Invalid UniProt alternative-products Sequence value in "
                    f"record {record.number}: {sequence_value}"
                )
            identifiers = [
                isoform_id.strip()
                for isoform_id in iso_ids.group(1).split(",")
                if isoform_id.strip()
            ]
            if identifiers:
                record.isoforms.append((identifiers, name.strip(), sequence_value))


def _validate_canonical_fasta(path: Path, index: sqlite3.Connection) -> None:
    index.execute("CREATE TEMP TABLE seen(accession TEXT PRIMARY KEY)")
    for identifier, sequence in _iter_fasta(path, role="canonical"):
        row = index.execute(
            "SELECT sequence FROM canonical WHERE accession=?", (identifier,)
        ).fetchone()
        if row is None or row[0] != sequence:
            raise ValueError(f"Canonical FASTA does not agree with DAT: {identifier}")
        try:
            index.execute("INSERT INTO seen VALUES (?)", (identifier,))
        except sqlite3.IntegrityError as error:
            raise ValueError(
                f"Duplicate canonical FASTA identifier: {identifier}"
            ) from error
    missing = index.execute(
        "SELECT accession FROM canonical EXCEPT SELECT accession FROM seen LIMIT 1"
    ).fetchone()
    if missing is not None:
        raise ValueError(f"Canonical FASTA is missing DAT accession: {missing[0]}")


def _write_isoform_fasta(
    path: Path,
    writer: csv.DictWriter[str],
    index: sqlite3.Connection,
) -> None:
    for fasta_identifier, sequence in _iter_fasta(path, role="isoform"):
        definitions = index.execute(
            "SELECT i.primary_accession, i.isoform_id "
            "FROM isoform_identifier ii "
            "JOIN isoform i USING (primary_accession, isoform_id) "
            "WHERE ii.identifier=? AND i.sequence_status='Alternative'",
            (fasta_identifier,),
        ).fetchall()
        if not definitions:
            raise ValueError(
                "Varsplic sequence has no materializable Alternative DAT "
                f"definition: {fasta_identifier}"
            )
        if len(definitions) != 1:
            raise ValueError(
                "Varsplic sequence has ambiguous Alternative owners: "
                f"{fasta_identifier}"
            )
        primary, isoform_id = definitions[0]
        try:
            index.execute(
                "INSERT INTO seen_isoform_sequence VALUES (?, ?)",
                (primary, isoform_id),
            )
        except sqlite3.IntegrityError as error:
            raise ValueError(
                f"Duplicate varsplic sequence product: {primary}/{isoform_id}"
            ) from error
        writer.writerow(
            {
                "sequence_id": isoform_id,
                "primary_accession": primary,
                "sequence_type": "isoform",
                "sequence": sequence,
                "length": len(sequence),
                "crc64": None,
                "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            }
        )
        index.execute(
            "UPDATE isoform SET sequence_id=? "
            "WHERE primary_accession=? AND isoform_id=? "
            "AND sequence_status='Alternative'",
            (isoform_id, primary, isoform_id),
        )


def _iter_fasta(path: Path, *, role: str) -> Iterator[tuple[str, str]]:
    with _open_fasta_text(path, role=role) as handle:
        identifier: str | None = None
        parts: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if identifier is not None:
                    yield (
                        identifier,
                        _validated_fasta_sequence(
                            parts, role=role, identifier=identifier
                        ),
                    )
                token = line[1:].split(None, 1)[0]
                header = token.split("|")
                if len(header) != 3 or header[0] != "sp":
                    raise ValueError(
                        f"UniProt {role} FASTA header must use sp|ID|ENTRY: {token}"
                    )
                identifier = header[1]
                if role == "isoform" and "-" not in identifier:
                    raise ValueError(
                        f"UniProt isoform FASTA identifier is not an IsoId: {identifier}"
                    )
                parts = []
            else:
                if identifier is None:
                    raise ValueError(
                        "Declared FASTA input does not begin with a header"
                    )
                sequence_line = line.rstrip("\r\n")
                if not sequence_line or any(
                    character.isspace() for character in sequence_line
                ):
                    raise ValueError(
                        f"UniProt {role} FASTA record has invalid sequence "
                        f"characters: {identifier}"
                    )
                parts.append(sequence_line)
        if identifier is not None:
            yield (
                identifier,
                _validated_fasta_sequence(parts, role=role, identifier=identifier),
            )
        else:
            raise ValueError("Declared FASTA input contains no records")


@contextmanager
def _open_fasta_text(path: Path, *, role: str) -> Generator[IO[str]]:
    try:
        with _open_text(path) as handle:
            yield handle
    except (gzip.BadGzipFile, EOFError, UnicodeDecodeError, zlib.error) as exc:
        if _media_type(path) == "application/gzip":
            raise ValueError(
                f"UniProt {role} FASTA input has an invalid gzip stream: {path}"
            ) from exc
        raise


def _validated_fasta_sequence(
    parts: Iterable[str], *, role: str, identifier: str
) -> str:
    sequence = "".join(parts)
    if not sequence:
        raise ValueError(
            f"UniProt {role} FASTA record has an empty sequence: {identifier}"
        )
    if any(residue not in _RESIDUE_MOLECULAR_WEIGHT for residue in sequence):
        raise ValueError(
            f"UniProt {role} FASTA record has invalid sequence characters: {identifier}"
        )
    return sequence


def _calculate_crc64(sequence: str) -> str:
    checksum = 0
    for residue in sequence:
        checksum = _CRC64_TABLE[(checksum ^ ord(residue)) & 0xFF] ^ (checksum >> 8)
    return f"{checksum:016X}"


def _calculate_molecular_weight(sequence: str, *, legacy_expasy: bool = False) -> int:
    weight = _MOLECULAR_WEIGHT_WATER
    for residue in sequence:
        if legacy_expasy and residue in _LEGACY_EXPASY_MOLECULAR_WEIGHT:
            weight += _LEGACY_EXPASY_MOLECULAR_WEIGHT[residue]
        else:
            weight += _RESIDUE_MOLECULAR_WEIGHT[residue]
    return (weight + 5000) // 10000


def _write_isoform_relations(
    writers: Mapping[str, csv.DictWriter[str]],
    index: sqlite3.Connection,
) -> None:
    missing = index.execute(
        "SELECT iv.variation_id FROM isoform_variation iv "
        "LEFT JOIN variation v USING (variation_id) "
        "WHERE v.variation_id IS NULL LIMIT 1"
    ).fetchone()
    if missing is not None:
        raise ValueError(f"Isoform references an absent VAR_SEQ feature: {missing[0]}")
    missing_sequence = index.execute(
        "SELECT primary_accession, isoform_id FROM isoform "
        "WHERE sequence_status='Alternative' AND sequence_id IS NULL LIMIT 1"
    ).fetchone()
    has_varsplic = (
        index.execute("SELECT count(*) FROM seen_isoform_sequence").fetchone() or (0,)
    )[0]
    if has_varsplic and missing_sequence is not None:
        raise ValueError(
            "Varsplic FASTA is missing Alternative isoform: "
            f"{missing_sequence[0]}/{missing_sequence[1]}"
        )
    for isoform_id, primary, name, isoform_order, status, sequence_id in index.execute(
        "SELECT isoform_id, primary_accession, name, isoform_order, "
        "sequence_status, sequence_id "
        "FROM isoform ORDER BY primary_accession, isoform_order"
    ):
        writers["protein_isoform"].writerow(
            {
                "isoform_id": isoform_id,
                "primary_accession": primary,
                "name": name,
                "isoform_order": isoform_order,
                "sequence_status": status,
                "sequence_id": sequence_id,
            }
        )
    for primary, isoform_id, identifier, identifier_order, is_main in index.execute(
        "SELECT primary_accession, isoform_id, identifier, identifier_order, is_main "
        "FROM isoform_identifier ORDER BY primary_accession, isoform_id, identifier_order"
    ):
        writers["protein_isoform_identifier"].writerow(
            {
                "primary_accession": primary,
                "isoform_id": isoform_id,
                "identifier": identifier,
                "identifier_order": identifier_order,
                "is_main": bool(is_main),
            }
        )
    for primary, isoform_id, variation_id, variation_order in index.execute(
        "SELECT primary_accession, isoform_id, variation_id, variation_order "
        "FROM isoform_variation "
        "ORDER BY primary_accession, isoform_id, variation_order"
    ):
        writers["protein_isoform_variation"].writerow(
            {
                "primary_accession": primary,
                "isoform_id": isoform_id,
                "variation_id": variation_id,
                "variation_order": variation_order,
            }
        )


@contextmanager
def _open_text(path: Path) -> Generator[IO[str]]:
    with path.open("rb") as raw:
        magic = raw.read(6)
    if magic.startswith(b"\x1f\x8b"):
        opener = gzip.open
    elif magic.startswith(b"BZh"):
        opener = bz2.open
    elif magic.startswith(b"\xfd7zXZ"):
        opener = lzma.open
    else:
        opener = Path.open
    handle = opener(path, "rt", encoding="utf-8", errors="strict")
    try:
        yield handle
    finally:
        handle.close()


def _media_type(path: Path) -> str:
    with path.open("rb") as handle:
        magic = handle.read(6)
    if magic.startswith(b"\x1f\x8b"):
        return "application/gzip"
    if magic.startswith(b"BZh"):
        return "application/x-bzip2"
    if magic.startswith(b"\xfd7zXZ"):
        return "application/x-xz"
    return "text/plain"
