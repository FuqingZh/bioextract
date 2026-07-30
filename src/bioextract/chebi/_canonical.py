from __future__ import annotations

import gzip
import io
import re
import tarfile
import zipfile
from collections import Counter
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import polars as pl

from bioextract._publication import RelationSpec, ValidationIssue

_CHEBI_ID = re.compile(r"^CHEBI:[0-9]+$")
_QUOTED = re.compile(r'"((?:[^"\\]|\\.)*)"')
_PROPERTY_NAME = re.compile(r"(?:[/#]|chebi/)([A-Za-z0-9_]+)$")
_RELATION_NAMES = {
    "BFO:0000051": "has_part",
    "RO:0000087": "has_role",
    "RO:0018033": "is_conjugate_base_of",
    "RO:0018034": "is_conjugate_acid_of",
    "RO:0018036": "is_tautomer_of",
    "RO:0018037": "is_substituent_group_from",
    "RO:0018038": "has_functional_parent",
    "RO:0018039": "is_enantiomer_of",
    "RO:0018040": "has_parent_hydride",
}


class ChEBIIntegrityError(RuntimeError):
    """Raised when canonical ChEBI entities cannot be published safely."""


@dataclass(frozen=True, slots=True)
class CanonicalBuild:
    relations: tuple[RelationSpec, ...]
    validation_issues: tuple[ValidationIssue, ...]


def build_canonical_relations(
    file_obo: Path,
    *,
    file_sdf: Path | None = None,
) -> CanonicalBuild:
    frames, issues = _read_chebi_obo(file_obo)
    if file_sdf is not None:
        structure, structure_issues = _read_sdf(
            file_sdf,
            canonical_ids=set(frames["compound"]["chebi_id"].to_list()),
        )
        frames["compound_structure"] = structure
        issues.extend(structure_issues)
    else:
        frames["compound_structure"] = _empty_structure_frame()

    roles = {
        "compound": "entity",
        "secondary_id": "identifier_mapping",
        "compound_name": "annotation",
        "compound_cross_reference": "cross_reference",
        "compound_relation": "relationship",
        "compound_structure": "structure",
        "compound_wurcs": "structure",
    }
    return CanonicalBuild(
        relations=tuple(
            RelationSpec(name, frame.lazy(), role=roles[name])
            for name, frame in frames.items()
        ),
        validation_issues=tuple(issues),
    )


def _read_chebi_obo(
    file_obo: Path,
) -> tuple[dict[str, pl.DataFrame], list[ValidationIssue]]:
    compounds: list[dict[str, object]] = []
    secondary_ids: list[dict[str, str]] = []
    names: list[dict[str, str]] = []
    xrefs: list[dict[str, str]] = []
    relations: list[dict[str, str]] = []
    wurcs: list[dict[str, str]] = []
    stanza: dict[str, list[str]] | None = None
    stanza_kind: str | None = None
    typedef_names: dict[str, str] = {}

    def flush() -> None:
        nonlocal stanza
        if stanza is None or "id" not in stanza:
            return
        stanza_id = stanza["id"][0].strip()
        if stanza_kind == "Typedef":
            name = _first(stanza, "name")
            if name:
                typedef_names[stanza_id] = _snake(name)
            return
        if stanza_kind != "Term":
            return
        if not _CHEBI_ID.fullmatch(stanza_id):
            raise ChEBIIntegrityError(
                f"Canonical ChEBI term has invalid id: {stanza_id!r}"
            )
        properties = _properties(stanza.get("property_value", []))
        subsets = stanza.get("subset", [])
        compounds.append(
            {
                "chebi_id": stanza_id,
                "preferred_name": _first(stanza, "name"),
                "definition": _quoted(_first(stanza, "def")),
                "star_rating": _star_rating(subsets),
                "is_obsolete": _first(stanza, "is_obsolete") == "true",
                "formula": _property(
                    properties, "formula", "generalized_empirical_formula"
                ),
                "charge": _optional_int(_property(properties, "charge")),
                "average_mass": _optional_float(
                    _property(properties, "mass", "average_mass")
                ),
                "monoisotopic_mass": _optional_float(
                    _property(properties, "monoisotopicmass", "monoisotopic_mass")
                ),
                "smiles": _property(properties, "smiles", "smiles_string"),
                "inchi": _property(properties, "inchi", "inchi_string"),
                "inchi_key": _property(
                    properties, "inchikey", "inchi_key", "inchi_key_string"
                ),
            }
        )
        for secondary_id in stanza.get("alt_id", []):
            secondary_ids.append(
                {
                    "secondary_chebi_id": secondary_id.strip(),
                    "chebi_id": stanza_id,
                }
            )
        for synonym in stanza.get("synonym", []):
            names.append(
                {
                    "chebi_id": stanza_id,
                    "name": _quoted(synonym) or "",
                    "scope": _synonym_scope(synonym),
                }
            )
        for xref in stanza.get("xref", []):
            xref_id = xref.split(maxsplit=1)[0].strip()
            prefix, separator, accession = xref_id.partition(":")
            if not separator:
                prefix, accession = xref_id, xref_id
            xrefs.append(
                {
                    "chebi_id": stanza_id,
                    "source_prefix": _xref_prefix(prefix),
                    "accession": accession.strip(),
                    "xref_id": xref_id,
                }
            )
        for parent in stanza.get("is_a", []):
            relations.append(
                {
                    "subject_chebi_id": stanza_id,
                    "relation_type": "is_a",
                    "relation_id": "is_a",
                    "object_chebi_id": parent.split()[0],
                }
            )
        for relationship in stanza.get("relationship", []):
            fields = relationship.split()
            if len(fields) >= 2:
                relation_id = fields[0]
                relations.append(
                    {
                        "subject_chebi_id": stanza_id,
                        "relation_type": typedef_names.get(
                            relation_id,
                            _RELATION_NAMES.get(relation_id, _snake(relation_id)),
                        ),
                        "relation_id": relation_id,
                        "object_chebi_id": fields[1],
                    }
                )
        for value in _property_values(properties, "wurcs", "wurcs_representation"):
            wurcs.append({"chebi_id": stanza_id, "wurcs": value})

    with _open_text(file_obo, preferred_suffix=".obo") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("[") and line.endswith("]"):
                flush()
                stanza_kind = line[1:-1]
                stanza = {}
                continue
            if stanza is None or not line or line.startswith("!"):
                continue
            key, separator, value = line.partition(":")
            if separator:
                stanza.setdefault(key, []).append(value.strip())
        flush()

    ids = [str(row["chebi_id"]) for row in compounds]
    if not ids:
        raise ChEBIIntegrityError("Core ChEBI OBO contains no canonical compounds")
    duplicates = sorted(value for value, count in Counter(ids).items() if count > 1)
    if duplicates:
        raise ChEBIIntegrityError(
            f"Core ChEBI OBO contains duplicate canonical IDs: {duplicates[:5]}"
        )
    canonical_ids = set(ids)
    issues: list[ValidationIssue] = []
    secondary_ids = _retain_valid_dependents(
        secondary_ids,
        target_field="chebi_id",
        canonical_ids=canonical_ids,
        source_name="chebi_obo",
        relation_name="secondary_id",
        identifier_namespace="chebi",
        identifier_field="secondary_chebi_id",
        issues=issues,
    )
    names = _retain_valid_dependents(
        names,
        target_field="chebi_id",
        canonical_ids=canonical_ids,
        source_name="chebi_obo",
        relation_name="compound_name",
        identifier_namespace="chebi",
        identifier_field="chebi_id",
        issues=issues,
    )
    xrefs = _retain_valid_dependents(
        xrefs,
        target_field="chebi_id",
        canonical_ids=canonical_ids,
        source_name="chebi_obo",
        relation_name="compound_cross_reference",
        identifier_namespace_field="source_prefix",
        identifier_field="accession",
        issues=issues,
    )
    wurcs = _retain_valid_dependents(
        wurcs,
        target_field="chebi_id",
        canonical_ids=canonical_ids,
        source_name="chebi_obo",
        relation_name="compound_wurcs",
        identifier_namespace="chebi",
        identifier_field="chebi_id",
        issues=issues,
    )
    relations = _retain_valid_relations(relations, canonical_ids, issues)

    return (
        {
            "compound": pl.DataFrame(
                compounds,
                schema={
                    "chebi_id": pl.String,
                    "preferred_name": pl.String,
                    "definition": pl.String,
                    "star_rating": pl.Int8,
                    "is_obsolete": pl.Boolean,
                    "formula": pl.String,
                    "charge": pl.Int32,
                    "average_mass": pl.Float64,
                    "monoisotopic_mass": pl.Float64,
                    "smiles": pl.String,
                    "inchi": pl.String,
                    "inchi_key": pl.String,
                },
            ),
            "secondary_id": pl.DataFrame(
                secondary_ids,
                schema={
                    "secondary_chebi_id": pl.String,
                    "chebi_id": pl.String,
                },
            ),
            "compound_name": pl.DataFrame(
                names,
                schema={
                    "chebi_id": pl.String,
                    "name": pl.String,
                    "scope": pl.String,
                },
            ),
            "compound_cross_reference": pl.DataFrame(
                xrefs,
                schema={
                    "chebi_id": pl.String,
                    "source_prefix": pl.String,
                    "accession": pl.String,
                    "xref_id": pl.String,
                },
            ),
            "compound_relation": pl.DataFrame(
                relations,
                schema={
                    "subject_chebi_id": pl.String,
                    "relation_type": pl.String,
                    "relation_id": pl.String,
                    "object_chebi_id": pl.String,
                },
            ),
            "compound_wurcs": pl.DataFrame(
                wurcs,
                schema={"chebi_id": pl.String, "wurcs": pl.String},
            ),
        },
        issues,
    )


def _read_sdf(
    file_sdf: Path,
    *,
    canonical_ids: set[str],
) -> tuple[pl.DataFrame, list[ValidationIssue]]:
    rows: list[dict[str, object]] = []
    issues: list[ValidationIssue] = []
    with _open_text(file_sdf, preferred_suffix=".sdf") as handle:
        for record_number, record in enumerate(_iter_sdf_records(handle), start=1):
            lines = record.splitlines()
            properties = _sdf_properties(lines)
            raw_id = properties.get("ChEBI ID") or properties.get("CHEBI ID")
            if raw_id is None:
                continue
            chebi_id = _canonical_chebi_id(raw_id)
            if chebi_id not in canonical_ids:
                issues.append(
                    _foreign_key_issue(
                        source_name="chebi_sdf",
                        relation_name="compound_structure",
                        identifier_namespace="chebi",
                        identifier_value=chebi_id,
                        referenced_identifier=chebi_id,
                        source_record_number=record_number,
                    )
                )
                continue
            marker = next(
                (index for index, line in enumerate(lines) if line.startswith("> <")),
                len(lines),
            )
            molfile = "\n".join(lines[:marker]).rstrip()
            rows.append(
                {
                    "chebi_id": chebi_id,
                    "structure_index": 1,
                    "molfile": molfile,
                }
            )
    return (
        pl.DataFrame(
            rows,
            schema={
                "chebi_id": pl.String,
                "structure_index": pl.Int32,
                "molfile": pl.String,
            },
        ),
        issues,
    )


def _retain_valid_dependents(
    rows: Sequence[dict[str, str]],
    *,
    target_field: str,
    canonical_ids: set[str],
    source_name: str,
    relation_name: str,
    identifier_field: str,
    issues: list[ValidationIssue],
    identifier_namespace: str | None = None,
    identifier_namespace_field: str | None = None,
) -> list[dict[str, str]]:
    retained: list[dict[str, str]] = []
    for record_number, row in enumerate(rows, start=1):
        target = row[target_field]
        if target in canonical_ids:
            retained.append(row)
            continue
        namespace = (
            row[identifier_namespace_field]
            if identifier_namespace_field is not None
            else identifier_namespace
        )
        issues.append(
            _foreign_key_issue(
                source_name=source_name,
                relation_name=relation_name,
                identifier_namespace=namespace,
                identifier_value=row.get(identifier_field),
                referenced_identifier=target,
                source_record_number=record_number,
            )
        )
    return retained


def _retain_valid_relations(
    rows: Sequence[dict[str, str]],
    canonical_ids: set[str],
    issues: list[ValidationIssue],
) -> list[dict[str, str]]:
    retained: list[dict[str, str]] = []
    for record_number, row in enumerate(rows, start=1):
        missing = next(
            (
                row[field]
                for field in ("subject_chebi_id", "object_chebi_id")
                if row[field] not in canonical_ids
            ),
            None,
        )
        if missing is None:
            retained.append(row)
            continue
        issues.append(
            _foreign_key_issue(
                source_name="chebi_obo",
                relation_name="compound_relation",
                identifier_namespace="chebi",
                identifier_value=row["subject_chebi_id"],
                referenced_identifier=missing,
                source_record_number=record_number,
            )
        )
    return retained


def _foreign_key_issue(
    *,
    source_name: str,
    relation_name: str,
    identifier_namespace: str | None,
    identifier_value: str | None,
    referenced_identifier: str | None,
    source_record_number: int,
) -> ValidationIssue:
    return ValidationIssue(
        severity="warning",
        issue_code="foreign_key_violation",
        source_name=source_name,
        relation_name=relation_name,
        identifier_namespace=identifier_namespace,
        identifier_value=identifier_value,
        referenced_relation="compound",
        referenced_identifier=referenced_identifier,
        source_record_number=source_record_number,
        message=(
            f"Skipped {relation_name} record referencing absent canonical "
            f"compound {referenced_identifier}"
        ),
    )


def _properties(values: Sequence[str]) -> dict[str, list[str]]:
    parsed: dict[str, list[str]] = {}
    for value in values:
        fields = value.split(maxsplit=1)
        if len(fields) != 2:
            continue
        raw_name, raw_value = fields
        match = _PROPERTY_NAME.search(raw_name)
        name = (match.group(1) if match else raw_name.rsplit(":", 1)[-1]).lower()
        quoted = _quoted(raw_value)
        parsed.setdefault(name, []).append(
            quoted if quoted is not None else raw_value.split()[0]
        )
    return parsed


def _property(
    properties: dict[str, list[str]],
    *names: str,
) -> str | None:
    values = _property_values(properties, *names)
    return values[0] if values else None


def _property_values(
    properties: dict[str, list[str]],
    *names: str,
) -> list[str]:
    for name in names:
        values = properties.get(name.lower())
        if values:
            return values
    return []


def _star_rating(subsets: Sequence[str]) -> int | None:
    for subset in subsets:
        match = re.search(r"([1-3])(?:_|:|\s*)STAR", subset, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def _xref_prefix(value: str) -> str:
    return re.sub(r"[ _]+", ".", value.strip().lower())


def _snake(value: str) -> str:
    normalized = re.sub(r"[^0-9A-Za-z]+", "_", value).strip("_").lower()
    return normalized or "related_to"


def _quoted(value: str | None) -> str | None:
    if value is None:
        return None
    match = _QUOTED.search(value)
    if match is None:
        return None
    return match.group(1).replace(r"\"", '"').replace(r"\\", "\\")


def _synonym_scope(value: str) -> str:
    quoted = _QUOTED.search(value)
    suffix = value[quoted.end() :] if quoted is not None else value
    match = re.search(r"\b(EXACT|BROAD|NARROW|RELATED)\b", suffix)
    return match.group(1) if match else "RELATED"


def _first(stanza: dict[str, list[str]], key: str) -> str | None:
    values = stanza.get(key)
    return values[0] if values else None


def _optional_int(value: str | None) -> int | None:
    try:
        return None if value is None else int(value)
    except ValueError:
        return None


def _optional_float(value: str | None) -> float | None:
    try:
        return None if value is None else float(value)
    except ValueError:
        return None


def _canonical_chebi_id(value: str) -> str:
    text = value.strip()
    if text.upper().startswith("CHEBI:"):
        text = text.split(":", maxsplit=1)[1]
    return f"CHEBI:{int(text)}"


def _iter_sdf_records(handle: TextIO) -> Iterator[str]:
    lines: list[str] = []
    for line in handle:
        if line.rstrip("\r\n") == "$$$$":
            yield "".join(lines)
            lines.clear()
        else:
            lines.append(line)
    if lines:
        yield "".join(lines)


def _sdf_properties(lines: Sequence[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    index = 0
    while index < len(lines):
        match = re.match(r">\s*<([^>]+)>", lines[index])
        if match is None:
            index += 1
            continue
        index += 1
        collected: list[str] = []
        while index < len(lines) and lines[index].strip():
            collected.append(lines[index].strip())
            index += 1
        values[match.group(1)] = "\n".join(collected)
    return values


def _empty_structure_frame() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "chebi_id": pl.String,
            "structure_index": pl.Int32,
            "molfile": pl.String,
        }
    )


@contextmanager
def _open_text(path: Path, *, preferred_suffix: str) -> Iterator[TextIO]:
    with path.open("rb") as probe:
        magic = probe.read(4)
    if magic.startswith(b"\x1f\x8b"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            names = [
                info.filename
                for info in archive.infolist()
                if not info.is_dir()
                and info.filename.lower().endswith(preferred_suffix)
            ]
            if len(names) != 1:
                raise ValueError(
                    f"Expected one {preferred_suffix} member in archive: {path}"
                )
            with (
                archive.open(names[0]) as raw,
                io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle,
            ):
                yield handle
        return
    if tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            members = [
                member
                for member in archive.getmembers()
                if member.isfile() and member.name.lower().endswith(preferred_suffix)
            ]
            if len(members) != 1:
                raise ValueError(
                    f"Expected one {preferred_suffix} member in archive: {path}"
                )
            raw = archive.extractfile(members[0])
            if raw is None:
                raise ValueError(f"Cannot read archive member: {members[0].name}")
            with (
                raw,
                io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle,
            ):
                yield handle
        return
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        yield handle
