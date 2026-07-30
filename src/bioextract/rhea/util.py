"""Streaming readers for Rhea RDF/XML, SDF, and release metadata."""

from __future__ import annotations

import gzip
import hashlib
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO
from xml.etree import ElementTree

_RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_RDFS_NS = "http://www.w3.org/2000/01/rdf-schema#"
_RHEA_NS = "http://rdf.rhea-db.org/"

_RDF_DESCRIPTION = f"{{{_RDF_NS}}}Description"
_RDF_ABOUT = f"{{{_RDF_NS}}}about"
_RDF_RESOURCE = f"{{{_RDF_NS}}}resource"
_SIDE_PATTERN = re.compile(r"^(?P<master_id>\d+)_(?P<side>[LR])$")


@dataclass(slots=True)
class RdfReaction:
    """One exact Rhea reaction record parsed from RDF/XML."""

    rhea_id: int
    accession: str | None = None
    equation: str | None = None
    equation_html: str | None = None
    status: str | None = None
    comment: str | None = None
    is_balanced: bool | None = None
    is_transport: bool | None = None
    side_ids: list[str] = field(default_factory=list)
    substrate_side_ids: list[str] = field(default_factory=list)
    product_side_ids: list[str] = field(default_factory=list)
    bidirectional_side_ids: list[str] = field(default_factory=list)
    citations: set[str] = field(default_factory=set)


@dataclass(slots=True)
class RdfSide:
    """One master-reaction side."""

    side_id: str
    master_id: int
    side: str
    curated_order: int | None = None


@dataclass(slots=True)
class RdfParticipant:
    """One RDF reaction participant."""

    participant_id: str
    compound_id: str | None = None
    location: str | None = None


@dataclass(slots=True)
class RdfCompound:
    """One Rhea RDF compound or reactive-part record."""

    compound_id: str
    rhea_compound_id: int | None = None
    source_accession: str | None = None
    name: str | None = None
    name_html: str | None = None
    formula: str | None = None
    charge: str | None = None
    chebi_id: int | None = None
    underlying_chebi_id: int | None = None
    polymerization_index: str | None = None
    position: str | None = None
    types: set[str] = field(default_factory=set)
    reactive_part_ids: set[str] = field(default_factory=set)


@dataclass(frozen=True, slots=True)
class RdfMembership:
    """One side-to-participant membership and coefficient predicate."""

    side_id: str
    participant_id: str
    coefficient_predicate: str


@dataclass(slots=True)
class RdfSnapshot:
    """Materialized semantic records from one streamed Rhea RDF/XML source."""

    reactions: dict[int, RdfReaction]
    sides: dict[str, RdfSide]
    participants: dict[str, RdfParticipant]
    compounds: dict[str, RdfCompound]
    memberships: set[RdfMembership]
    coefficient_by_predicate: dict[str, str]


@dataclass(frozen=True, slots=True)
class SdfStructure:
    """One structure record from ``rhea.sdf``."""

    accession: str
    role: str | None
    chebi_xref: str | None
    generic_compound_accession: str | None
    underlying_chebi_polymer_accession: str | None
    formula: str | None
    charge: str | None
    name: str | None
    molblock: str


def is_gzip_file(file_path: Path) -> bool:
    """Return whether a file starts with the gzip magic bytes."""

    with file_path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text_auto(file_path: Path) -> TextIO:
    """Open plain text or gzip-compressed text using content detection."""

    if is_gzip_file(file_path):
        return gzip.open(file_path, "rt", encoding="utf-8")
    return file_path.open("rt", encoding="utf-8")


def calculate_sha256(file_path: Path) -> str:
    """Calculate one source-file SHA-256 digest."""

    with file_path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def infer_media_type(file_path: Path) -> str:
    """Infer a compact source media type from its logical suffix."""

    name = file_path.name.removesuffix(".gz").lower()
    if name.endswith(".rdf"):
        return "application/rdf+xml"
    if name.endswith(".sdf"):
        return "chemical/x-mdl-sdfile"
    if name.endswith(".tsv"):
        return "text/tab-separated-values"
    if name.endswith(".properties"):
        return "text/x-java-properties"
    return "text/plain"


def read_release_properties(file_path: Path) -> tuple[int | None, str | None]:
    """Read Rhea release number and date from ``rhea-release.properties``."""

    values: dict[str, str] = {}
    with open_text_auto(file_path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", maxsplit=1)
            values[key.strip()] = value.strip()

    number_text = values.get("rhea.release.number")
    return (
        int(number_text) if number_text is not None else None,
        values.get("rhea.release.date"),
    )


def parse_rdf_snapshot(file_path: Path) -> RdfSnapshot:
    """Stream one Rhea RDF/XML file into normalized semantic records."""

    reactions: dict[int, RdfReaction] = {}
    sides: dict[str, RdfSide] = {}
    participants: dict[str, RdfParticipant] = {}
    compounds: dict[str, RdfCompound] = {}
    memberships_raw: set[RdfMembership] = set()
    coefficient_by_predicate: dict[str, str] = {}

    with open_text_auto(file_path) as handle:
        for _, element in ElementTree.iterparse(handle, events=("end",)):
            if element.tag != _RDF_DESCRIPTION:
                continue

            subject_uri = element.attrib.get(_RDF_ABOUT)
            if subject_uri is None or not subject_uri.startswith(_RHEA_NS):
                element.clear()
                continue

            subject = subject_uri.removeprefix(_RHEA_NS)
            if subject.isdigit():
                _update_reaction(reactions, int(subject), element)
            elif match_side := _SIDE_PATTERN.fullmatch(subject):
                _update_side_and_memberships(
                    sides,
                    memberships_raw,
                    subject,
                    int(match_side.group("master_id")),
                    match_side.group("side"),
                    element,
                )
            elif subject.startswith("Participant_"):
                _update_participant(participants, subject, element)
            elif subject.startswith("Compound_"):
                _update_compound(compounds, subject, element)
            elif subject.startswith("contains"):
                coefficient = _child_text(element, _RHEA_NS, "coefficient")
                if coefficient is not None:
                    coefficient_by_predicate[subject] = coefficient

            element.clear()

    memberships = _deduplicate_memberships(
        memberships_raw,
        coefficient_by_predicate=coefficient_by_predicate,
    )
    return RdfSnapshot(
        reactions=reactions,
        sides=sides,
        participants=participants,
        compounds=compounds,
        memberships=memberships,
        coefficient_by_predicate=coefficient_by_predicate,
    )


def iter_sdf_structures(file_path: Path) -> Iterator[SdfStructure]:
    """Yield structure records from plain or gzip-compressed Rhea SDF."""

    record_lines: list[str] = []
    with open_text_auto(file_path) as handle:
        for line in handle:
            if line.rstrip("\r\n") == "$$$$":
                if record_lines:
                    yield _parse_sdf_record(record_lines)
                    record_lines.clear()
                continue
            record_lines.append(line)

    if record_lines:
        yield _parse_sdf_record(record_lines)


def compound_type(compound: RdfCompound) -> str | None:
    """Resolve the most specific Rhea compound type."""

    precedence = (
        "ReactivePart",
        "Polymer",
        "GenericPolypeptide",
        "GenericPolynucleotide",
        "GenericHeteropolysaccharide",
        "GenericSmallMolecule",
        "GenericCompound",
        "SmallMolecule",
        "Compound",
    )
    return next((value for value in precedence if value in compound.types), None)


def public_compound_accession(source_accession: str | None) -> str | None:
    """Map internal GENERIC/POLYMER accessions to public RHEA-COMP accessions."""

    if source_accession is None:
        return None
    prefix, separator, identifier = source_accession.partition(":")
    if separator and prefix in {"GENERIC", "POLYMER"}:
        return f"RHEA-COMP:{identifier}"
    return source_accession


def reaction_direction(
    reaction: RdfReaction,
) -> tuple[int | None, str | None]:
    """Derive master ID and direction from official RDF side references.

    Obsolete RDF tombstones have no side references and therefore retain null
    direction metadata rather than receiving an invented direction.
    """

    if reaction.side_ids:
        return reaction.rhea_id, "UN"
    if reaction.substrate_side_ids:
        master_id, side = _split_side_id(reaction.substrate_side_ids[0])
        return master_id, "LR" if side == "L" else "RL"
    if reaction.bidirectional_side_ids:
        master_id, _ = _split_side_id(reaction.bidirectional_side_ids[0])
        return master_id, "BI"
    if (reaction.status or "").lower() == "obsolete":
        return None, None
    raise ValueError(f"Cannot derive direction for RHEA:{reaction.rhea_id}")


def coefficient_numeric(coefficient: str) -> float | None:
    """Convert a numeric coefficient while preserving symbolic values as null."""

    try:
        return float(coefficient)
    except ValueError:
        return None


def _update_reaction(
    reactions: dict[int, RdfReaction],
    rhea_id: int,
    element: ElementTree.Element,
) -> None:
    reaction = reactions.setdefault(rhea_id, RdfReaction(rhea_id=rhea_id))
    for child in element:
        namespace, name = _split_tag(child.tag)
        resource = child.attrib.get(_RDF_RESOURCE)
        text = child.text
        if namespace == _RHEA_NS:
            if name == "accession":
                reaction.accession = text
            elif name == "equation":
                reaction.equation = text
            elif name == "htmlEquation":
                reaction.equation_html = text
            elif name == "status" and resource is not None:
                reaction.status = _uri_tail(resource)
            elif name == "isChemicallyBalanced":
                reaction.is_balanced = _parse_bool(text)
            elif name == "isTransport":
                reaction.is_transport = _parse_bool(text)
            elif name == "side" and resource is not None:
                reaction.side_ids.append(_uri_tail(resource))
            elif name == "substrates" and resource is not None:
                reaction.substrate_side_ids.append(_uri_tail(resource))
            elif name == "products" and resource is not None:
                reaction.product_side_ids.append(_uri_tail(resource))
            elif name == "substratesOrProducts" and resource is not None:
                reaction.bidirectional_side_ids.append(_uri_tail(resource))
            elif name == "citation" and resource is not None:
                reaction.citations.add(_uri_tail(resource))
        elif namespace == _RDFS_NS and name == "comment":
            reaction.comment = text


def _update_side_and_memberships(
    sides: dict[str, RdfSide],
    memberships: set[RdfMembership],
    side_id: str,
    master_id: int,
    side: str,
    element: ElementTree.Element,
) -> None:
    record = sides.setdefault(
        side_id,
        RdfSide(side_id=side_id, master_id=master_id, side=side),
    )
    for child in element:
        namespace, name = _split_tag(child.tag)
        if namespace != _RHEA_NS:
            continue
        if name == "curatedOrder" and child.text is not None:
            record.curated_order = int(child.text)
        elif name.startswith("contains"):
            resource = child.attrib.get(_RDF_RESOURCE)
            if resource is not None:
                memberships.add(
                    RdfMembership(
                        side_id=side_id,
                        participant_id=_uri_tail(resource),
                        coefficient_predicate=name,
                    )
                )


def _update_participant(
    participants: dict[str, RdfParticipant],
    participant_id: str,
    element: ElementTree.Element,
) -> None:
    participant = participants.setdefault(
        participant_id,
        RdfParticipant(participant_id=participant_id),
    )
    for child in element:
        namespace, name = _split_tag(child.tag)
        resource = child.attrib.get(_RDF_RESOURCE)
        if namespace != _RHEA_NS or resource is None:
            continue
        if name == "compound":
            participant.compound_id = _uri_tail(resource)
        elif name == "location":
            participant.location = _uri_tail(resource).lower()


def _update_compound(
    compounds: dict[str, RdfCompound],
    compound_id: str,
    element: ElementTree.Element,
) -> None:
    compound = compounds.setdefault(
        compound_id,
        RdfCompound(compound_id=compound_id),
    )
    for child in element:
        namespace, name = _split_tag(child.tag)
        resource = child.attrib.get(_RDF_RESOURCE)
        text = child.text
        if namespace == _RHEA_NS:
            if name == "id" and text is not None:
                compound.rhea_compound_id = int(text)
            elif name == "accession":
                compound.source_accession = text
            elif name == "name":
                compound.name = text
            elif name == "htmlName":
                compound.name_html = text
            elif name == "formula":
                compound.formula = text
            elif name == "charge":
                compound.charge = text
            elif name == "chebi" and resource is not None:
                compound.chebi_id = _parse_chebi_uri(resource)
            elif name == "underlyingChebi" and resource is not None:
                compound.underlying_chebi_id = _parse_chebi_uri(resource)
            elif name == "polymerizationIndex":
                compound.polymerization_index = text
            elif name == "position":
                compound.position = text
            elif name == "reactivePart" and resource is not None:
                compound.reactive_part_ids.add(_uri_tail(resource))
        elif (
            namespace == _RDFS_NS
            and name == "subClassOf"
            and resource is not None
            and resource.startswith(_RHEA_NS)
        ):
            compound.types.add(_uri_tail(resource))


def _deduplicate_memberships(
    memberships: set[RdfMembership],
    *,
    coefficient_by_predicate: dict[str, str],
) -> set[RdfMembership]:
    membership_by_key: dict[tuple[str, str], RdfMembership] = {}
    for membership in sorted(
        memberships,
        key=lambda value: (
            value.side_id,
            value.participant_id,
            value.coefficient_predicate == "contains",
        ),
    ):
        key = (membership.side_id, membership.participant_id)
        existing = membership_by_key.get(key)
        if existing is None:
            membership_by_key[key] = membership
            continue
        if (
            existing.coefficient_predicate == "contains"
            and membership.coefficient_predicate in coefficient_by_predicate
        ):
            membership_by_key[key] = membership
    return set(membership_by_key.values())


def _parse_sdf_record(lines: list[str]) -> SdfStructure:
    properties: dict[str, str] = {}
    index = 0
    while index < len(lines):
        line = lines[index].rstrip("\r\n")
        if line.startswith("> <") and line.endswith(">"):
            key = line[3:-1]
            values: list[str] = []
            index += 1
            while index < len(lines) and lines[index].strip():
                values.append(lines[index].rstrip("\r\n"))
                index += 1
            properties[key] = "\n".join(values)
        index += 1

    end_molblock = next(
        (index for index, line in enumerate(lines) if line.rstrip("\r\n") == "M  END"),
        None,
    )
    if end_molblock is None:
        raise ValueError("Rhea SDF record has no 'M  END' marker")

    accession = properties.get("ACCESSION") or lines[0].strip()
    charge_text = properties.get("Charge")
    return SdfStructure(
        accession=accession,
        role=properties.get("ROLE"),
        chebi_xref=properties.get("CHEBI_XREF"),
        generic_compound_accession=properties.get("GENERIC_COMPOUND"),
        underlying_chebi_polymer_accession=properties.get("UNDERLYING_CHEBI_POLYMER"),
        formula=properties.get("Formula"),
        charge=charge_text,
        name=properties.get("Rhea_ascii_name"),
        molblock="".join(lines[: end_molblock + 1]).rstrip() + "\n",
    )


def _child_text(
    element: ElementTree.Element,
    namespace: str,
    name: str,
) -> str | None:
    child = element.find(f"{{{namespace}}}{name}")
    return None if child is None else child.text


def _split_tag(tag: str) -> tuple[str, str]:
    if not tag.startswith("{"):
        return "", tag
    namespace, name = tag[1:].split("}", maxsplit=1)
    return namespace, name


def _uri_tail(uri: str) -> str:
    return uri.rsplit("/", maxsplit=1)[-1]


def _parse_chebi_uri(uri: str) -> int:
    return int(_uri_tail(uri).removeprefix("CHEBI_"))


def _parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    return value.strip().lower() == "true"


def _split_side_id(side_id: str) -> tuple[int, str]:
    match = _SIDE_PATTERN.fullmatch(side_id)
    if match is None:
        raise ValueError(f"Invalid Rhea side identifier: {side_id}")
    return int(match.group("master_id")), match.group("side")
