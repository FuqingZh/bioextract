from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NoReturn, cast

import polars as pl
import pyarrow as _pyarrow  # pyright: ignore[reportMissingTypeStubs]

from bioextract._publication import SourceFileRecord, ValidationIssue
from bioextract.errors import IntegrityError

from .constant import (
    MEDIA_TYPE_TSV,
    SCHEMA_GENE_ANNOTATION,
    SCHEMA_KO_ANNOTATION,
    SCHEMA_ORGANISM,
)

_KO_ID = re.compile(r"^K[0-9]{5}$")
_ORGANISM_PATHWAY = re.compile(r"^([a-z]{3,4})([0-9]{5})$")
_GLOBAL_PATHWAY = re.compile(r"^(map|ko)([0-9]{5})$")
_GENOME_ID = re.compile(r"^T[0-9]+$")
pa: Any = _pyarrow


@dataclass(frozen=True, slots=True)
class ParsedFile:
    rows: tuple[tuple[str, ...], ...]
    source: SourceFileRecord
    issues: tuple[ValidationIssue, ...]


@dataclass(frozen=True, slots=True)
class OrganismAggregate:
    organism: pl.DataFrame
    genes: pl.DataFrame
    ko_ids: frozenset[str]
    sources: tuple[SourceFileRecord, ...]
    issues: tuple[ValidationIssue, ...]


def parse_organism_list(
    path: Path | None,
) -> tuple[
    dict[str, dict[str, object]],
    tuple[SourceFileRecord, ...],
    tuple[ValidationIssue, ...],
]:
    if path is None:
        return {}, (), ()
    parsed = _read_rows(
        path,
        logical_name="global/organism_list",
        columns=4,
    )
    result: dict[str, dict[str, object]] = {}
    for line_number, (genome_id, code, name, taxonomy) in enumerate(parsed.rows, 1):
        if _GENOME_ID.fullmatch(genome_id) is None:
            _raise_row(path, "global/organism_list", line_number, "invalid genome_id")
        value: dict[str, object] = {
            "genome_id": genome_id,
            "organism_name": name,
            "taxonomy_lineage": [
                item for raw in taxonomy.split(";") if (item := raw.strip())
            ],
        }
        previous = result.get(code)
        if previous is not None and previous != value:
            _raise_row(
                path,
                "global/organism_list",
                line_number,
                f"conflicting metadata for organism {code!r}",
            )
        result[code] = value
    return result, (parsed.source,), parsed.issues


def aggregate_organism(
    organism_code: str,
    roles: Mapping[str, Path],
    *,
    organism_metadata: Mapping[str, object] | None,
    organism_list_available: bool,
) -> OrganismAggregate:
    parsed: dict[str, ParsedFile] = {}
    specifications = {
        "gene_list": 4,
        "uniprot_conversion": 2,
        "ncbi_gene_conversion": 2,
        "gene_ko": 2,
        "gene_pathway": 2,
    }
    for role, path in roles.items():
        columns = specifications[role]
        parsed[role] = _read_rows(
            path,
            logical_name=f"organism/{organism_code}/{role}",
            columns=columns,
        )

    genes: dict[str, dict[str, object]] = {}

    def gene(gene_id: str, *, role: str, line_number: int) -> dict[str, object]:
        _validate_gene_id(
            gene_id,
            organism_code=organism_code,
            path=roles[role],
            role=role,
            line_number=line_number,
        )
        return genes.setdefault(gene_id, {})

    for line_number, row in enumerate(parsed.get("gene_list", _EMPTY).rows, 1):
        gene_id, gene_type, genomic_position, display = row
        target = gene(gene_id, role="gene_list", line_number=line_number)
        attributes = _parse_gene_display(display)
        attributes.update(
            gene_type=gene_type or None,
            genomic_position=genomic_position or None,
        )
        previous = target.get("attributes")
        if previous is not None and previous != attributes:
            _raise_row(
                roles["gene_list"],
                f"organism/{organism_code}/gene_list",
                line_number,
                f"conflicting attributes for gene {gene_id!r}",
            )
        target["attributes"] = attributes

    for line_number, (xref, gene_id) in enumerate(
        parsed.get("uniprot_conversion", _EMPTY).rows, 1
    ):
        prefix, separator, accession = xref.partition(":")
        if not separator or not prefix or not accession:
            _raise_row(
                roles["uniprot_conversion"],
                f"organism/{organism_code}/uniprot_conversion",
                line_number,
                "invalid UniProt conversion identifier",
            )
        target = gene(gene_id, role="uniprot_conversion", line_number=line_number)
        cast("set[str]", target.setdefault("uniprot", set[str]())).add(accession)

    for line_number, (xref, gene_id) in enumerate(
        parsed.get("ncbi_gene_conversion", _EMPTY).rows, 1
    ):
        prefix, separator, ncbi_id = xref.partition(":")
        if prefix != "ncbi-geneid" or not separator or not ncbi_id.isdigit():
            _raise_row(
                roles["ncbi_gene_conversion"],
                f"organism/{organism_code}/ncbi_gene_conversion",
                line_number,
                "invalid NCBI Gene conversion identifier",
            )
        target = gene(gene_id, role="ncbi_gene_conversion", line_number=line_number)
        cast("set[str]", target.setdefault("ncbi_gene", set[str]())).add(ncbi_id)

    ko_ids: set[str] = set()
    for line_number, (gene_id, raw_ko) in enumerate(
        parsed.get("gene_ko", _EMPTY).rows, 1
    ):
        ko_id = _strip_required_prefix(
            raw_ko,
            "ko:",
            path=roles["gene_ko"],
            role=f"organism/{organism_code}/gene_ko",
            line_number=line_number,
        )
        if _KO_ID.fullmatch(ko_id) is None:
            _raise_row(
                roles["gene_ko"],
                f"organism/{organism_code}/gene_ko",
                line_number,
                "invalid KO identifier",
            )
        target = gene(gene_id, role="gene_ko", line_number=line_number)
        cast("set[str]", target.setdefault("ko", set[str]())).add(ko_id)
        ko_ids.add(ko_id)

    for line_number, (gene_id, raw_pathway) in enumerate(
        parsed.get("gene_pathway", _EMPTY).rows, 1
    ):
        pathway_id = _strip_required_prefix(
            raw_pathway,
            "path:",
            path=roles["gene_pathway"],
            role=f"organism/{organism_code}/gene_pathway",
            line_number=line_number,
        )
        match = _ORGANISM_PATHWAY.fullmatch(pathway_id)
        if match is None or match.group(1) != organism_code:
            _raise_row(
                roles["gene_pathway"],
                f"organism/{organism_code}/gene_pathway",
                line_number,
                "pathway identifier does not match organism",
            )
        target = gene(gene_id, role="gene_pathway", line_number=line_number)
        cast(
            "set[tuple[str, str]]",
            target.setdefault("pathway", set[tuple[str, str]]()),
        ).add((pathway_id, f"map{match.group(2)}"))

    gene_columns: dict[str, list[object]] = {
        name: [] for name in SCHEMA_GENE_ANNOTATION
    }
    for gene_id in sorted(genes):
        values = genes[gene_id]
        attributes = values.get("attributes")
        attrs: Mapping[str, object] = (
            cast("Mapping[str, object]", attributes)
            if isinstance(attributes, dict)
            else {}
        )
        output_row: dict[str, object] = {
            "organism_code": organism_code,
            "kegg_gene_id": gene_id,
            "gene_type": attrs.get("gene_type"),
            "genomic_position": attrs.get("genomic_position"),
            "gene_symbol": attrs.get("gene_symbol"),
            "gene_aliases": _optional_scalar_list(
                values,
                attrs,
                key="gene_aliases",
                capability="gene_list" in roles,
            ),
            "gene_description": attrs.get("gene_description"),
            "uniprot_mappings": _struct_list(
                values.get("uniprot"),
                field="uniprot_id",
                capability="uniprot_conversion" in roles,
            ),
            "ncbi_gene_mappings": _struct_list(
                values.get("ncbi_gene"),
                field="ncbi_gene_id",
                capability="ncbi_gene_conversion" in roles,
            ),
            "ko_mappings": _struct_list(
                values.get("ko"),
                field="ko_id",
                capability="gene_ko" in roles,
            ),
            "pathway_mappings": (
                None
                if "gene_pathway" not in roles
                else [
                    (pathway_id, map_id)
                    for pathway_id, map_id in sorted(
                        cast(
                            "set[tuple[str, str]]",
                            values.get("pathway", set[tuple[str, str]]()),
                        )
                    )
                ]
            ),
        }
        for name in SCHEMA_GENE_ANNOTATION:
            gene_columns[name].append(output_row[name])
    metadata = organism_metadata or {}
    organism_row: dict[str, object] = {
        "organism_code": organism_code,
        "genome_id": metadata.get("genome_id"),
        "organism_name": metadata.get("organism_name"),
        "taxonomy_lineage": (
            metadata.get("taxonomy_lineage")
            if metadata
            else ([] if organism_list_available else None)
        ),
    }
    return OrganismAggregate(
        organism=_dataframe_from_columns(
            {name: [organism_row[name]] for name in SCHEMA_ORGANISM},
            SCHEMA_ORGANISM,
        ),
        genes=_dataframe_from_columns(gene_columns, SCHEMA_GENE_ANNOTATION),
        ko_ids=frozenset(ko_ids),
        sources=tuple(value.source for value in parsed.values()),
        issues=tuple(issue for value in parsed.values() for issue in value.issues),
    )


def aggregate_kos(
    ko_ids: set[str],
    path: Path | None,
) -> tuple[pl.DataFrame, tuple[SourceFileRecord, ...], tuple[ValidationIssue, ...]]:
    pathways: dict[str, set[tuple[str, str, str]]] = {}
    sources: tuple[SourceFileRecord, ...] = ()
    issues: tuple[ValidationIssue, ...] = ()
    if path is not None:
        parsed = _read_rows(
            path,
            logical_name="global/ko_pathway",
            columns=2,
        )
        sources = (parsed.source,)
        issues = parsed.issues
        for line_number, (raw_ko, raw_pathway) in enumerate(parsed.rows, 1):
            ko_id = _strip_required_prefix(
                raw_ko,
                "ko:",
                path=path,
                role="global/ko_pathway",
                line_number=line_number,
            )
            pathway_id = _strip_required_prefix(
                raw_pathway,
                "path:",
                path=path,
                role="global/ko_pathway",
                line_number=line_number,
            )
            match = _GLOBAL_PATHWAY.fullmatch(pathway_id)
            if _KO_ID.fullmatch(ko_id) is None or match is None:
                _raise_row(
                    path,
                    "global/ko_pathway",
                    line_number,
                    "invalid KO-to-pathway edge",
                )
            pathways.setdefault(ko_id, set()).add(
                (pathway_id, match.group(1), f"map{match.group(2)}")
            )
            ko_ids.add(ko_id)
    ko_columns: dict[str, list[object]] = {name: [] for name in SCHEMA_KO_ANNOTATION}
    for ko_id in sorted(ko_ids):
        ko_columns["ko_id"].append(ko_id)
        ko_columns["pathway_mappings"].append(
            None
            if path is None
            else [
                (pathway_id, namespace, map_id)
                for pathway_id, namespace, map_id in sorted(pathways.get(ko_id, set()))
            ]
        )
    return (
        _dataframe_from_columns(ko_columns, SCHEMA_KO_ANNOTATION),
        sources,
        issues,
    )


def _dataframe_from_columns(
    columns: Mapping[str, Sequence[object]], schema: Mapping[str, object]
) -> pl.DataFrame:
    """Build nested frames through Arrow without per-cell dict conversion."""
    arrays = {
        name: pa.array(values, type=_arrow_type(schema[name]))
        for name, values in columns.items()
    }
    return cast(
        "pl.DataFrame",
        pl.from_arrow(pa.table(arrays), rechunk=False),  # pyright: ignore[reportUnknownMemberType]
    )


def _arrow_type(dtype: object) -> Any:
    if dtype == pl.String:
        return pa.string()
    if isinstance(dtype, pl.List):
        return pa.list_(_arrow_type(dtype.inner))
    if isinstance(dtype, pl.Struct):
        return pa.struct(
            [(field.name, _arrow_type(field.dtype)) for field in dtype.fields]
        )
    raise AssertionError(f"Unsupported KEGG mapping dtype: {dtype}")


def _read_rows(
    path: Path,
    *,
    logical_name: str,
    columns: int,
) -> ParsedFile:
    try:
        content = path.read_bytes()
    except OSError as error:
        raise IntegrityError(f"Cannot read {logical_name}: {path}") from error
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as error:
        raise IntegrityError(f"Invalid UTF-8 in {logical_name}: {path}") from error
    rows: list[tuple[str, ...]] = []
    seen: set[tuple[str, ...]] = set()
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        fields = tuple(field.strip() for field in line.split("\t"))
        if len(fields) != columns:
            _raise_row(
                path,
                logical_name,
                line_number,
                f"expected {columns} columns, found {len(fields)}",
            )
        if fields in seen:
            continue
        seen.add(fields)
        rows.append(fields)
    source = SourceFileRecord(
        logical_name=logical_name,
        path=path,
        media_type=MEDIA_TYPE_TSV,
        bytes=len(content),
        sha256=None,
    )
    return ParsedFile(tuple(rows), source, ())


def _parse_gene_display(display: str) -> dict[str, object]:
    if ";" not in display:
        return {
            "gene_symbol": None,
            "gene_aliases": [],
            "gene_description": display.strip() or None,
        }
    names, description = display.split(";", 1)
    parsed_names = [name for raw in names.split(",") if (name := raw.strip())]
    return {
        "gene_symbol": parsed_names[0] if parsed_names else None,
        "gene_aliases": sorted(set(parsed_names[1:])),
        "gene_description": description.strip() or None,
    }


def _optional_scalar_list(
    values: Mapping[str, object],
    attributes: Mapping[str, object],
    *,
    key: str,
    capability: bool,
) -> object:
    del values
    return attributes.get(key, []) if capability else None


def _struct_list(value: object, *, field: str, capability: bool) -> object:
    if not capability:
        return None
    items = cast("set[str]", value) if isinstance(value, set) else set[str]()
    return [(item,) for item in sorted(items)]


def _strip_required_prefix(
    value: str,
    prefix: str,
    *,
    path: Path,
    role: str,
    line_number: int,
) -> str:
    if not value.startswith(prefix) or len(value) == len(prefix):
        _raise_row(path, role, line_number, f"expected prefix {prefix!r}")
    return value[len(prefix) :]


def _validate_gene_id(
    gene_id: str,
    *,
    organism_code: str,
    path: Path,
    role: str,
    line_number: int,
) -> None:
    if (
        not gene_id.startswith(f"{organism_code}:")
        or len(gene_id) <= len(organism_code) + 1
    ):
        _raise_row(
            path,
            f"organism/{organism_code}/{role}",
            line_number,
            f"gene identifier does not match organism {organism_code!r}",
        )


def _raise_row(path: Path, role: str, line_number: int, message: str) -> NoReturn:
    raise IntegrityError(
        f"KEGG mapping parse error: role={role!r}, path={path}, "
        f"line={line_number}: {message}"
    )


_EMPTY = ParsedFile(
    rows=(),
    source=SourceFileRecord("empty", Path("."), MEDIA_TYPE_TSV, 0),
    issues=(),
)
