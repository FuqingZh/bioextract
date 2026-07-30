from __future__ import annotations

import hashlib
import json
import re
import tarfile
import tempfile
import zipfile
from collections.abc import Mapping
from contextlib import ExitStack
from pathlib import Path
from typing import Any, Literal

import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    ValidationIssue,
    write_duckdb_publication,
)

from .core import (
    ENTRY_ROLES,
    RELATION_ROLES,
    SCHEMA_VERSION,
    MetabolicSnapshot,
    ModuleDefinitionParser,
    ast_rows,
    db_links,
    discover_release_layout,
    entry_id,
    iter_records,
    joined,
    numeric,
    open_text,
    parse_equation,
    record_names,
)

_RELATION_COLUMNS = {
    "compound_pubchem": ("compound_id", "pubchem_id"),
    "compound_reaction": ("compound_id", "reaction_id"),
    "reaction_enzyme": ("reaction_id", "ec_number"),
    "reaction_ko": ("reaction_id", "ko_id"),
    "reaction_module": ("reaction_id", "module_id"),
    "reaction_pathway": ("reaction_id", "pathway_id"),
    "module_pathway": ("module_id", "pathway_id"),
}
_EC_NUMBER = re.compile(r"(?<![\d.])\d+\.\d+\.\d+\.[\d-]+(?![\d.])")


def _safe_target(root: Path, name: str) -> Path:
    target = (root / name).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError(f"Unsafe path in KEGG release archive: {name}")
    return target


def extract_archive(path: Path, root: Path) -> None:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for member in archive.infolist():
                _safe_target(root, member.filename)
                if member.is_dir():
                    continue
                target = _safe_target(root, member.filename)
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(member) as source, target.open("wb") as output:
                    while chunk := source.read(1024 * 1024):
                        output.write(chunk)
        return
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            _safe_target(root, member.name)
            if member.issym() or member.islnk():
                raise ValueError(
                    f"Links are not allowed in KEGG release archives: {member.name}"
                )
            if not member.isfile():
                continue
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(f"Cannot read KEGG archive member: {member.name}")
            target = _safe_target(root, member.name)
            target.parent.mkdir(parents=True, exist_ok=True)
            with source, target.open("wb") as output:
                while chunk := source.read(1024 * 1024):
                    output.write(chunk)


class _Spool:
    def __init__(self, root: Path):
        self.root = root
        self.stack = ExitStack()
        self.handles: dict[str, Any] = {}
        self.counts: dict[str, int] = {}

    def add(self, relation: str, row: Mapping[str, Any]) -> None:
        handle = self.handles.get(relation)
        if handle is None:
            path = self.root / f"{relation}.ndjson"
            handle = self.stack.enter_context(path.open("w", encoding="utf-8"))
            self.handles[relation] = handle
        handle.write(json.dumps(dict(row), separators=(",", ":")) + "\n")
        self.counts[relation] = self.counts.get(relation, 0) + 1

    def close(self) -> None:
        self.stack.close()

    def relations(self) -> list[RelationSpec]:
        return [
            RelationSpec(
                name,
                pl.scan_ndjson(self.root / f"{name}.ndjson")
                .unique()
                .sort(
                    pl.scan_ndjson(self.root / f"{name}.ndjson")
                    .collect_schema()
                    .names()
                ),
                role="canonical",
            )
            for name in sorted(self.counts)
        ]


def _read_list_ids(paths: tuple[Path, ...]) -> set[str]:
    result: set[str] = set()
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                if line.strip():
                    result.add(line.split("\t", 1)[0].split(":", 1)[-1].strip())
    return result


def _enzyme_replacement_targets(
    record: Mapping[str, list[str]],
    ec_number: str,
) -> tuple[str, ...]:
    entry = joined(record, "ENTRY") or ""
    definition = joined(record, "DEFINITION") or ""
    is_obsolete = "Obsolete" in entry.split()
    has_legacy_transfer = definition == "Transferred entry" or bool(
        record.get("TRANSFER")
    )
    if not (is_obsolete or has_legacy_transfer):
        return ()
    text = " ".join(
        (
            joined(record, "NAME") or "",
            joined(record, "COMMENT") or "",
            joined(record, "TRANSFER") or "",
        )
    )
    return tuple(sorted(set(_EC_NUMBER.findall(text)) - {ec_number}))


def _spool_entries(
    snapshot: MetabolicSnapshot, spool: _Spool
) -> tuple[
    dict[str, set[str]],
    set[tuple[str, str]],
    list[ValidationIssue],
]:
    ids: dict[str, set[str]] = {family: set() for family in ENTRY_ROLES}
    equation_pairs: set[tuple[str, str]] = set()
    issues: list[ValidationIssue] = []
    enzyme_replacements: list[tuple[str, str]] = []
    for family in ENTRY_ROLES:
        for record in iter_records(snapshot.sources.get(f"{family}_entries", ())):
            ident = entry_id(record)
            if ident in ids[family]:
                raise ValueError(f"Duplicate KEGG {family} primary ID: {ident}")
            ids[family].add(ident)
            names = record_names(record)
            if family == "compound":
                spool.add(
                    "compound",
                    {
                        "compound_id": ident,
                        "name": names[0] if names else None,
                        "formula": joined(record, "FORMULA"),
                        "exact_mass": numeric(joined(record, "EXACT_MASS")),
                        "molecular_weight": numeric(joined(record, "MOL_WEIGHT")),
                    },
                )
                for position, name in enumerate(names, 1):
                    spool.add(
                        "compound_name",
                        {
                            "compound_id": ident,
                            "position": position,
                            "name": name,
                            "is_primary": position == 1,
                        },
                    )
                for namespace, xref in db_links(record):
                    external = (
                        f"CHEBI:{xref.removeprefix('CHEBI:')}"
                        if namespace.casefold() == "chebi"
                        else xref
                    )
                    spool.add(
                        "compound_cross_reference",
                        {
                            "compound_id": ident,
                            "namespace": namespace.casefold(),
                            "external_id": external,
                            "relationship": "cross_reference",
                        },
                    )
            elif family == "reaction":
                equation = joined(record, "EQUATION")
                spool.add(
                    "reaction",
                    {
                        "reaction_id": ident,
                        "name": names[0] if names else None,
                        "definition": joined(record, "DEFINITION"),
                        "equation": equation,
                        "is_reversible": bool(equation and "<=>" in equation),
                    },
                )
                for position, name in enumerate(names, 1):
                    spool.add(
                        "reaction_name",
                        {
                            "reaction_id": ident,
                            "position": position,
                            "name": name,
                            "is_primary": position == 1,
                        },
                    )
                if equation:
                    for row in parse_equation(ident, equation):
                        if (
                            row["participant_namespace"] == "kegg_compound"
                            and "compound_entries" in snapshot.sources
                            and row["participant_id"] not in ids["compound"]
                        ):
                            issues.append(
                                ValidationIssue(
                                    "warning",
                                    "foreign_key_violation",
                                    "reaction_entries",
                                    "reaction_participant",
                                    "kegg_reaction",
                                    ident,
                                    "compound",
                                    row["participant_id"],
                                    None,
                                    "Reaction participant compound is absent",
                                )
                            )
                            continue
                        spool.add("reaction_participant", row)
                        if row["participant_namespace"] == "kegg_compound":
                            equation_pairs.add((row["participant_id"], ident))
                for namespace, xref in db_links(record):
                    external = (
                        f"RHEA:{xref.removeprefix('RHEA:')}"
                        if namespace.casefold() == "rhea"
                        else xref
                    )
                    spool.add(
                        "reaction_cross_reference",
                        {
                            "reaction_id": ident,
                            "namespace": namespace.casefold(),
                            "external_id": external,
                            "relationship": "cross_reference",
                        },
                    )
                for line in record.get("RCLASS", ()):
                    values = line.split(maxsplit=1)
                    spool.add(
                        "reaction_class",
                        {
                            "reaction_id": ident,
                            "rclass_id": values[0],
                            "compound_pair": values[1] if len(values) > 1 else None,
                        },
                    )
            elif family == "enzyme":
                definition = joined(record, "DEFINITION")
                entry = joined(record, "ENTRY") or ""
                replacements = _enzyme_replacement_targets(record, ident)
                is_obsolete = "Obsolete" in entry.split()
                status = (
                    "transferred"
                    if replacements
                    else "deleted"
                    if is_obsolete or definition == "Deleted entry"
                    else "active"
                )
                spool.add(
                    "enzyme",
                    {
                        "ec_number": ident,
                        "status": status,
                        "class": joined(record, "CLASS"),
                        "systematic_name": joined(record, "SYSNAME"),
                        "comment": joined(record, "COMMENT"),
                        "history": joined(record, "HISTORY"),
                    },
                )
                for position, name in enumerate(names, 1):
                    spool.add(
                        "enzyme_name",
                        {
                            "ec_number": ident,
                            "position": position,
                            "name": name,
                            "is_primary": position == 1,
                        },
                    )
                for line in record.get("ORTHOLOGY", ()):
                    for ko_id in re.findall(r"K\d{5}", line):
                        spool.add("enzyme_ko", {"ec_number": ident, "ko_id": ko_id})
                enzyme_replacements.extend(
                    (ident, replacement) for replacement in replacements
                )
                for namespace, xref in db_links(record):
                    spool.add(
                        "enzyme_cross_reference",
                        {
                            "ec_number": ident,
                            "namespace": namespace.casefold(),
                            "external_id": xref,
                            "relationship": "cross_reference",
                        },
                    )
            else:
                definition = joined(record, "DEFINITION")
                entry_parts = (joined(record, "ENTRY") or "").split()
                diagram = joined(record, "DIAGRAM")
                spool.add(
                    "module",
                    {
                        "module_id": ident,
                        "type": entry_parts[1] if len(entry_parts) > 1 else None,
                        "name": joined(record, "NAME"),
                        "class": joined(record, "CLASS"),
                        "definition": definition,
                        "diagram_id": diagram.split()[0] if diagram else None,
                    },
                )
                if definition:
                    for row in ast_rows(
                        ident, ModuleDefinitionParser(definition).parse()
                    ):
                        spool.add("module_definition_node", row)
                    for ko_id in sorted(set(re.findall(r"K\d{5}", definition))):
                        spool.add("module_ko", {"module_id": ident, "ko_id": ko_id})
                step = 0
                for line in record.get("REACTION", ()):
                    values = line.split(maxsplit=1)
                    if not values:
                        continue
                    for reaction_id in re.findall(r"R\d{5}", values[0]):
                        step += 1
                        spool.add(
                            "module_reaction_step",
                            {
                                "module_id": ident,
                                "position": step,
                                "reaction_id": reaction_id,
                                "equation": values[1] if len(values) > 1 else None,
                            },
                        )
                for line in record.get("COMPOUND", ()):
                    for compound_id in re.findall(r"C\d{5}", line):
                        spool.add(
                            "module_compound",
                            {"module_id": ident, "compound_id": compound_id},
                        )
    for ec_number, replacement in enzyme_replacements:
        if replacement not in ids["enzyme"]:
            issues.append(
                ValidationIssue(
                    "warning",
                    "foreign_key_violation",
                    "enzyme_entries",
                    "enzyme_replacement",
                    "ec",
                    ec_number,
                    "enzyme",
                    replacement,
                    None,
                    "Replacement EC number is absent",
                )
            )
            continue
        spool.add(
            "enzyme_replacement",
            {
                "ec_number": ec_number,
                "replacement_ec_number": replacement,
            },
        )
    return ids, equation_pairs, issues


def _spool_relations(
    snapshot: MetabolicSnapshot,
    spool: _Spool,
    ids: Mapping[str, set[str]],
    equation_pairs: set[tuple[str, str]],
) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    official_pairs: set[tuple[str, str]] = set()
    endpoint = {
        "compound_pubchem": (("compound_id", "compound"),),
        "compound_reaction": (("compound_id", "compound"), ("reaction_id", "reaction")),
        "reaction_enzyme": (("reaction_id", "reaction"), ("ec_number", "enzyme")),
        "reaction_ko": (("reaction_id", "reaction"),),
        "reaction_module": (("reaction_id", "reaction"), ("module_id", "module")),
        "reaction_pathway": (("reaction_id", "reaction"),),
        "module_pathway": (("module_id", "module"),),
    }
    for role in RELATION_ROLES:
        columns = _RELATION_COLUMNS[role]
        for path in snapshot.sources.get(role, ()):
            with open_text(path) as handle:
                for record_number, line in enumerate(handle, 1):
                    if not line.strip():
                        continue
                    values = [
                        value.split(":", 1)[-1].strip()
                        for value in line.rstrip().split("\t")[:2]
                    ]
                    row = dict(zip(columns, values, strict=True))
                    invalid = next(
                        (
                            (column, family)
                            for column, family in endpoint[role]
                            if f"{family}_entries" in snapshot.sources
                            and row[column] not in ids[family]
                        ),
                        None,
                    )
                    if invalid:
                        column, family = invalid
                        issues.append(
                            ValidationIssue(
                                "warning",
                                "foreign_key_violation",
                                role,
                                role,
                                column,
                                row[column],
                                family,
                                row[column],
                                record_number,
                                "Dependent relation endpoint is absent",
                            )
                        )
                        continue
                    spool.add(role, row)
                    if role == "compound_reaction":
                        official_pairs.add((row["compound_id"], row["reaction_id"]))
    if (
        "reaction_entries" in snapshot.sources
        and official_pairs
        and official_pairs != equation_pairs
    ):
        for compound_id, reaction_id in sorted(official_pairs ^ equation_pairs):
            issues.append(
                ValidationIssue(
                    "warning",
                    "equation_relation_mismatch",
                    "compound_reaction",
                    "reaction_participant",
                    "kegg_compound",
                    compound_id,
                    "reaction",
                    reaction_id,
                    None,
                    "Equation C-number participants disagree with compound_reaction.tsv",
                )
            )
    return issues


def _resolve_entry_archives(
    snapshot: MetabolicSnapshot, root: Path
) -> MetabolicSnapshot:
    sources = {role: list(paths) for role, paths in snapshot.sources.items()}
    for family in ENTRY_ROLES:
        role = f"{family}_entries"
        resolved: list[Path] = []
        for index, path in enumerate(sources.get(role, [])):
            if not (zipfile.is_zipfile(path) or tarfile.is_tarfile(path)):
                resolved.append(path)
                continue
            target = root / f"{role}-{index:04d}"
            target.mkdir()
            extract_archive(path, target)
            batches = sorted(target.rglob("*.keg"))
            if not batches:
                raise ValueError(
                    f"KEGG {role} archive contains no .keg entry batches: {path}"
                )
            resolved.extend(batches)
        if resolved:
            sources[role] = resolved
    return MetabolicSnapshot(
        {role: tuple(paths) for role, paths in sources.items()},
        snapshot.release_version,
        snapshot.complete_release,
    )


def _source_media_type(path: Path) -> str:
    if zipfile.is_zipfile(path):
        return "application/zip"
    if tarfile.is_tarfile(path):
        return "application/x-tar"
    if path.suffix == ".gz":
        return "application/gzip"
    return "text/plain"


def _source_records(
    snapshot: MetabolicSnapshot,
    *,
    include_source_hashes: bool,
) -> list[SourceFileRecord]:
    logical_sources: Mapping[str, tuple[Path, ...]]
    if snapshot.archive is not None:
        logical_sources = {"release_archive": (snapshot.archive,)}
    else:
        logical_sources = snapshot.sources
    records: list[SourceFileRecord] = []
    for role, paths in logical_sources.items():
        for index, source in enumerate(paths):
            digest = (
                hashlib.sha256(source.read_bytes()).hexdigest()
                if include_source_hashes
                else None
            )
            records.append(
                SourceFileRecord(
                    f"{role}:{index}",
                    source,
                    _source_media_type(source),
                    digest,
                )
            )
    return records


def publish(
    snapshot: MetabolicSnapshot,
    path: Path,
    *,
    if_exists: Literal["fail", "replace"],
    include_source_hashes: bool,
) -> DuckDBWriteResult:
    with tempfile.TemporaryDirectory(prefix="bioextract-kegg-metabolic-") as directory:
        root = Path(directory)
        effective = snapshot
        if snapshot.archive is not None:
            release = root / "release"
            release.mkdir()
            extract_archive(snapshot.archive, release)
            discovered_sources = discover_release_layout(release)
            missing = [
                role
                for role in (
                    *[f"{family}_entries" for family in ENTRY_ROLES],
                    *RELATION_ROLES,
                )
                if role not in discovered_sources
            ]
            if missing:
                raise ValueError(
                    "Incomplete KEGG metabolic release archive; "
                    f"missing roles: {missing}"
                )
            effective = MetabolicSnapshot(
                discovered_sources, snapshot.release_version, True
            )
        effective = _resolve_entry_archives(effective, root)
        spool = _Spool(root / "relations")
        spool.root.mkdir()
        ids, equation_pairs, issues = _spool_entries(effective, spool)
        if effective.complete_release:
            for family in ENTRY_ROLES:
                listed = _read_list_ids(effective.sources.get(f"{family}_list", ()))
                if listed != ids[family]:
                    raise ValueError(f"KEGG {family} list/entry mismatch")
        issues.extend(_spool_relations(effective, spool, ids, equation_pairs))
        spool.close()
        sources = _source_records(
            snapshot,
            include_source_hashes=include_source_hashes,
        )
        return write_duckdb_publication(
            spool.relations(),
            path,
            resource_name="kegg",
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="kegg-metabolic-flat-files-v1",
            sources=sources,
            scope="metabolic",
            release_version=effective.release_version,
            if_exists=if_exists,
            validation_issues=issues,
            extra_metadata={
                "bioextract.capabilities": ",".join(
                    sorted(
                        role
                        for role in effective.sources
                        if role.endswith("_entries") or role in RELATION_ROLES
                    )
                )
            },
        )
