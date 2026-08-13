from __future__ import annotations

import copy
import gzip
import os
import re
import tarfile
import zipfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Literal

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import register_deferred_frame_source
from bioextract._publication import (
    DuckDBWriteResult,
    validate_duckdb_metadata_v1,
)
from bioextract._shared import validate_group_ids
from bioextract.errors import CapabilityError

SCHEMA_VERSION = "kegg-metabolic-v0.1"
SOURCE_SCHEMA_PROFILE = "kegg-metabolic-flat-files-v1"
METADATA_VERSION = "2"
NAMESPACES = (
    "kegg_compound",
    "chebi",
    "pubchem",
    "kegg_reaction",
    "rhea",
    "ec",
    "ko",
    "kegg_module",
    "kegg_pathway",
)
type KEGGMetabolicNamespace = Literal[
    "kegg_compound",
    "chebi",
    "pubchem",
    "kegg_reaction",
    "rhea",
    "ec",
    "ko",
    "kegg_module",
    "kegg_pathway",
]
ENTRY_ROLES = ("compound", "reaction", "enzyme", "module")
RELATION_ROLES = (
    "compound_pubchem",
    "compound_reaction",
    "reaction_enzyme",
    "reaction_ko",
    "reaction_module",
    "reaction_pathway",
    "module_pathway",
)
LIST_ROLES = tuple(f"{family}_list" for family in ENTRY_ROLES)
ENTRY_SOURCE_ROLES = tuple(f"{family}_entries" for family in ENTRY_ROLES)
REQUIRED_SOURCE_ROLES = (*ENTRY_SOURCE_ROLES, *RELATION_ROLES)
SOURCE_ROLES = (*LIST_ROLES, *ENTRY_SOURCE_ROLES, *RELATION_ROLES)
_RELATION_ROLE_BY_FILENAME = {f"{role}.tsv": role for role in RELATION_ROLES}
_CAPABILITY_TABLES: Mapping[str, frozenset[str]] = {
    "compound_entries": frozenset(
        {"compound", "compound_name", "compound_cross_reference"}
    ),
    "reaction_entries": frozenset(
        {
            "reaction",
            "reaction_name",
            "reaction_participant",
            "reaction_cross_reference",
            "reaction_class",
        }
    ),
    "enzyme_entries": frozenset(
        {
            "enzyme",
            "enzyme_name",
            "enzyme_cross_reference",
            "enzyme_replacement",
            "enzyme_ko",
        }
    ),
    "module_entries": frozenset(
        {
            "module",
            "module_definition_node",
            "module_ko",
            "module_reaction_step",
            "module_compound",
        }
    ),
    **{role: frozenset({role}) for role in RELATION_ROLES},
}


@dataclass(frozen=True, slots=True)
class MetabolicSnapshot:
    """Resolved local inputs for one KEGG metabolic publication build.

    Examples:
        >>> snapshot = MetabolicSnapshot(sources={}, release_version="2026-07")
        >>> snapshot.release_version
        '2026-07'
    """

    sources: Mapping[str, tuple[Path, ...]]
    release_version: str | None = None
    complete_release: bool = False
    archive: Path | None = None


@dataclass(frozen=True, slots=True)
class MetabolicPublication:
    """Validated identity and inventory of one KEGG metabolic DuckDB.

    Examples:
        >>> publication = MetabolicPublication(
        ...     Path("kegg.duckdb"), {}, frozenset(), frozenset()
        ... )
        >>> publication.path.name
        'kegg.duckdb'
    """

    path: Path
    metadata: Mapping[str, str]
    tables: frozenset[str]
    capabilities: frozenset[str]


def _paths(
    value: os.PathLike[str] | str | Sequence[os.PathLike[str] | str] | None,
    *,
    role: str,
) -> tuple[Path, ...]:
    if value is None:
        return ()
    values = (value,) if isinstance(value, (str, os.PathLike)) else tuple(value)
    result: list[Path] = []
    for raw in values:
        path = Path(raw)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.is_dir():
            pattern = "*.keg" if role.endswith("_entries") else "*"
            result.extend(sorted(p for p in path.rglob(pattern) if p.is_file()))
        else:
            result.append(path)
    paths = tuple(sorted(result))
    if not paths:
        raise ValueError(
            f"KEGG metabolic input role {role!r} must resolve to at least one file"
        )
    return paths


_METABOLIC_TABLE_SCHEMAS: Mapping[str, SchemaDict] = {
    "reaction": {
        "reaction_id": pl.String,
        "name": pl.String,
        "definition": pl.String,
        "equation": pl.String,
        "is_reversible": pl.Boolean,
    },
    "reaction_participant": {
        "reaction_id": pl.String,
        "side": pl.String,
        "position": pl.Int64,
        "participant_namespace": pl.String,
        "participant_id": pl.String,
        "coefficient_text": pl.String,
        "coefficient_numeric": pl.Float64,
    },
    "reaction_enzyme": {"reaction_id": pl.String, "ec_number": pl.String},
    "reaction_ko": {"reaction_id": pl.String, "ko_id": pl.String},
    "reaction_module": {"reaction_id": pl.String, "module_id": pl.String},
    "reaction_pathway": {"reaction_id": pl.String, "pathway_id": pl.String},
    "compound": {
        "compound_id": pl.String,
        "name": pl.String,
        "formula": pl.String,
        "exact_mass": pl.Float64,
        "molecular_weight": pl.Float64,
    },
    "compound_cross_reference": {
        "compound_id": pl.String,
        "namespace": pl.String,
        "external_id": pl.String,
        "relationship": pl.String,
    },
    "reaction_cross_reference": {
        "reaction_id": pl.String,
        "namespace": pl.String,
        "external_id": pl.String,
        "relationship": pl.String,
    },
}


def _table_columns(publication: MetabolicPublication, table: str) -> list[str]:
    """Return the contract columns for a generated metabolic relation."""
    del publication
    try:
        return list(_METABOLIC_TABLE_SCHEMAS[table])
    except KeyError as error:
        raise KeyError(f"Unknown KEGG metabolic table: {table}") from error


def _selection_schema(columns: Iterable[str]) -> SchemaDict:
    types: dict[str, Any] = {
        "group_id": pl.String,
        "input_id": pl.String,
        "input_namespace": pl.String,
        "anchor_type": pl.String,
        "anchor_id": pl.String,
        "match_type": pl.String,
        "entity_type": pl.String,
        "entity_id": pl.String,
        "reason": pl.String,
    }
    for table_schema in _METABOLIC_TABLE_SCHEMAS.values():
        for name, dtype in table_schema.items():
            if name not in types:
                types[name] = dtype
    return {column: types.get(column, pl.String) for column in columns}


def from_metabolic_files(
    source: os.PathLike[str] | str | None = None,
    *,
    compound_list: os.PathLike[str] | str | None = None,
    compound_entries: os.PathLike[str]
    | str
    | Sequence[os.PathLike[str] | str]
    | None = None,
    reaction_list: os.PathLike[str] | str | None = None,
    reaction_entries: os.PathLike[str]
    | str
    | Sequence[os.PathLike[str] | str]
    | None = None,
    enzyme_list: os.PathLike[str] | str | None = None,
    enzyme_entries: os.PathLike[str]
    | str
    | Sequence[os.PathLike[str] | str]
    | None = None,
    module_list: os.PathLike[str] | str | None = None,
    module_entries: os.PathLike[str]
    | str
    | Sequence[os.PathLike[str] | str]
    | None = None,
    compound_pubchem: os.PathLike[str] | str | None = None,
    compound_reaction: os.PathLike[str] | str | None = None,
    reaction_enzyme: os.PathLike[str] | str | None = None,
    reaction_ko: os.PathLike[str] | str | None = None,
    reaction_module: os.PathLike[str] | str | None = None,
    reaction_pathway: os.PathLike[str] | str | None = None,
    module_pathway: os.PathLike[str] | str | None = None,
    release_version: str | None = None,
) -> MetabolicSnapshot:
    """Resolve explicit roles and one optional release source into a snapshot."""
    values: dict[str, Any] = locals()
    explicit_sources = {
        role: _paths(values[role], role=role)
        for role in SOURCE_ROLES
        if values[role] is not None
    }
    validate_source_inventory(explicit_sources)

    if source is None:
        if explicit_sources:
            return MetabolicSnapshot(
                sources=explicit_sources,
                release_version=release_version,
            )
        raise ValueError("At least one KEGG metabolic input must be provided")

    source_path = Path(source)
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    if source_path.is_dir():
        discovered_sources = discover_release_layout(source_path)
        final_sources = dict(discovered_sources)
        final_sources.update(explicit_sources)
        validate_complete_release(final_sources)
        validate_source_inventory(final_sources)
        return MetabolicSnapshot(
            sources=final_sources,
            release_version=release_version,
            complete_release=True,
        )

    archive_roles = _inspect_release_archive(source_path)
    validate_source_inventory({"release_archive": (source_path,), **explicit_sources})
    missing = sorted(
        role
        for role in REQUIRED_SOURCE_ROLES
        if role not in archive_roles and role not in explicit_sources
    )
    if missing:
        raise ValueError(
            f"Incomplete KEGG metabolic release archive; missing roles: {missing}"
        )
    return MetabolicSnapshot(
        sources=explicit_sources,
        release_version=release_version,
        complete_release=True,
        archive=source_path,
    )


def discover_release_layout(root: Path) -> dict[str, tuple[Path, ...]]:
    """Resolve exactly one KEGG raw layout without traversal-order precedence."""
    candidates: list[Path] = []
    candidate_ids: set[Path] = set()
    for candidate in (root, *sorted(root.rglob("raw"))):
        if not candidate.is_dir():
            continue
        identity = candidate.resolve()
        if identity in candidate_ids:
            continue
        candidate_ids.add(identity)
        candidates.append(candidate)

    layouts = [
        (candidate, sources)
        for candidate in candidates
        if (sources := _discover_release_candidate(candidate)) is not None
    ]
    if not layouts:
        raise ValueError(f"KEGG metabolic source contains no release layout: {root}")
    if len(layouts) > 1:
        raise ValueError(
            "KEGG metabolic source contains multiple release layouts: "
            + ", ".join(str(candidate) for candidate, _ in layouts)
        )
    return layouts[0][1]


def _discover_release_candidate(
    raw: Path,
) -> dict[str, tuple[Path, ...]] | None:
    found: dict[str, tuple[Path, ...]] = {}
    recognized = False
    for family in ENTRY_ROLES:
        file_list = raw / family / "list.tsv"
        dir_entries = raw / family / "entries"
        if file_list.is_file():
            recognized = True
            found[f"{family}_list"] = (file_list,)
        if dir_entries.is_dir():
            recognized = True
            entries = tuple(sorted(dir_entries.glob("*.keg")))
            if entries:
                found[f"{family}_entries"] = entries
    for role in RELATION_ROLES:
        relation = raw / f"{role}.tsv"
        if relation.is_file():
            recognized = True
            found[role] = (relation,)
    return found if recognized else None


def validate_complete_release(
    sources: Mapping[str, tuple[Path, ...]],
) -> None:
    """Reject a source-backed inventory missing a required metabolic role."""
    missing = sorted(role for role in REQUIRED_SOURCE_ROLES if not sources.get(role))
    if missing:
        raise ValueError(f"Incomplete KEGG metabolic release; missing roles: {missing}")


def validate_source_inventory(
    sources: Mapping[str, tuple[Path, ...]],
) -> None:
    """Reject empty roles and physical-file reuse in a final inventory."""
    physical_roles: dict[tuple[int, int], tuple[str, Path]] = {}
    for role, paths in sources.items():
        if not paths:
            raise ValueError(
                f"KEGG metabolic input role {role!r} must resolve to at least one file"
            )
        for path in paths:
            resolved = path.resolve()
            stat = resolved.stat()
            physical_id = (stat.st_dev, stat.st_ino)
            if physical_id in physical_roles:
                other_role, other_path = physical_roles[physical_id]
                raise ValueError(
                    f"KEGG metabolic input roles {other_role!r} and {role!r} "
                    f"refer to the same physical file: {other_path}"
                )
            physical_roles[physical_id] = (role, resolved)


def _inspect_release_archive(path: Path) -> frozenset[str]:
    members: list[str]
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members = [
                member.filename for member in archive.infolist() if not member.is_dir()
            ]
    elif tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            links = [
                member.name
                for member in archive.getmembers()
                if member.issym() or member.islnk()
            ]
            if links:
                raise ValueError(
                    "Links are not allowed in KEGG release archives: "
                    + ", ".join(links)
                )
            members = [
                member.name for member in archive.getmembers() if member.isfile()
            ]
    else:
        raise ValueError(f"Unsupported KEGG metabolic release archive: {path}")

    layouts: dict[PurePosixPath, set[str]] = {}
    normalized_members = [PurePosixPath(name) for name in members]
    if len(normalized_members) != len(set(normalized_members)):
        raise ValueError(
            f"KEGG metabolic release archive has duplicate members: {path}"
        )
    for name, member in zip(members, normalized_members, strict=True):
        if member.is_absolute() or ".." in member.parts:
            raise ValueError(f"Unsafe path in KEGG release archive: {name}")
        match = _archive_member_role(member)
        if match is None:
            continue
        candidate, role = match
        layouts.setdefault(candidate, set()).add(role)

    if not layouts:
        raise ValueError(
            f"KEGG metabolic release archive contains no release layout: {path}"
        )
    if len(layouts) > 1:
        raise ValueError(
            "KEGG metabolic release archive contains multiple release layouts: "
            + ", ".join(str(candidate) for candidate in sorted(layouts, key=str))
        )
    return frozenset(next(iter(layouts.values())))


def _archive_member_role(
    member: PurePosixPath,
) -> tuple[PurePosixPath, str] | None:
    parts = member.parts
    candidate: PurePosixPath
    role: str
    if member.name in _RELATION_ROLE_BY_FILENAME:
        candidate = member.parent
        role = _RELATION_ROLE_BY_FILENAME[member.name]
    elif len(parts) >= 2 and parts[-2] in ENTRY_ROLES and member.name == "list.tsv":
        candidate = member.parent.parent
        role = f"{parts[-2]}_list"
    elif (
        len(parts) >= 3
        and parts[-3] in ENTRY_ROLES
        and parts[-2] == "entries"
        and member.suffix == ".keg"
    ):
        candidate = member.parent.parent.parent
        role = f"{parts[-3]}_entries"
    else:
        return None

    if candidate != PurePosixPath(".") and candidate.name != "raw":
        return None
    return candidate, role


def open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8", errors="replace")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8", errors="replace")
    )


def iter_records(paths: Iterable[Path]) -> Iterator[dict[str, list[str]]]:
    for path in paths:
        with open_text(path) as handle:
            record: dict[str, list[str]] = {}
            current: str | None = None
            for raw in handle:
                line = raw.rstrip("\r\n")
                if line == "///":
                    if record:
                        yield record
                    record, current = {}, None
                    continue
                label = line[:12].strip()
                value = line[12:].strip()
                if label:
                    current = label
                    record.setdefault(label, []).append(value)
                elif current is not None:
                    record[current].append(value)
            if record:
                raise ValueError(f"Unterminated KEGG record in {path}")


def joined(record: Mapping[str, list[str]], key: str, sep: str = " ") -> str | None:
    value = sep.join(record.get(key, ())).strip()
    return value or None


def record_names(record: Mapping[str, list[str]]) -> list[str]:
    raw = joined(record, "NAME") or ""
    return [x.strip() for x in raw.split(";") if x.strip()]


def db_links(record: Mapping[str, list[str]]) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    current: str | None = None
    for line in record.get("DBLINKS", ()):
        match = re.match(r"([^:]+):\s*(.*)", line)
        if match:
            current, values = match.groups()
        elif current:
            values = line
        else:
            continue
        for value in values.split():
            result.append(((current or "").strip(), value.strip()))
    return result


def entry_id(record: Mapping[str, list[str]]) -> str:
    entry = joined(record, "ENTRY")
    if not entry:
        raise ValueError("KEGG record is missing ENTRY")
    parts = entry.split()
    if parts[0] == "EC" and len(parts) > 1:
        return parts[1]
    return parts[0]


def numeric(value: str | None) -> float | None:
    try:
        return None if value is None else float(value)
    except ValueError:
        return None


def parse_equation(reaction_id: str, equation: str) -> list[dict[str, Any]]:
    if "<=>" not in equation:
        raise ValueError(f"Unparseable KEGG equation for {reaction_id}: {equation}")
    rows: list[dict[str, Any]] = []
    for side, text in zip(("left", "right"), equation.split("<=>", 1), strict=True):
        for position, item in enumerate(re.split(r"\s+\+\s+", text.strip()), 1):
            match = re.fullmatch(
                r"(?:(\([^)]*\)|[^\s]+)\s+)?([CG]\d{5})(?:\([^)]*\))?",
                item.strip(),
            )
            if not match:
                raise ValueError(f"Unparseable participant in {reaction_id}: {item}")
            coefficient, participant = match.groups()
            coefficient = coefficient or "1"
            rows.append(
                {
                    "reaction_id": reaction_id,
                    "side": side,
                    "position": position,
                    "participant_namespace": "kegg_compound"
                    if participant.startswith("C")
                    else "kegg_glycan",
                    "participant_id": participant,
                    "coefficient_text": coefficient,
                    "coefficient_numeric": numeric(coefficient),
                }
            )
    return rows


@dataclass(slots=True)
class _Node:
    kind: str
    children: list[_Node] = field(default_factory=lambda: [])
    value: str | None = None


def _tokenize(expression: str) -> list[str]:
    return re.findall(r"K\d{5}|M\d{5}|[()+,\-]|[^\s()+,\-]+", expression)


class ModuleDefinitionParser:
    def __init__(self, expression: str):
        self.tokens = _tokenize(expression)
        self.pos = 0

    def parse(self) -> _Node:
        node = self._sequence()
        if self.pos != len(self.tokens):
            raise ValueError(f"Unexpected module token: {self.tokens[self.pos]}")
        return node

    def _sequence(self) -> _Node:
        children = [self._alternative()]
        while self.pos < len(self.tokens) and self.tokens[self.pos] not in {")", ","}:
            children.append(self._alternative())
        return children[0] if len(children) == 1 else _Node("sequence", children)

    def _alternative(self) -> _Node:
        children = [self._complex()]
        while self._take(","):
            children.append(self._complex())
        return children[0] if len(children) == 1 else _Node("alternative", children)

    def _complex(self) -> _Node:
        children = [self._atom()]
        while self._take("+"):
            children.append(self._atom())
        return children[0] if len(children) == 1 else _Node("complex", children)

    def _atom(self) -> _Node:
        optional = self._take("-")
        if optional and self._take("-"):
            return _Node("optional")
        if self._take("("):
            node = self._sequence()
            if not self._take(")"):
                raise ValueError("Unclosed module definition group")
        elif self.pos < len(self.tokens) and re.fullmatch(
            r"[KM]\d{5}", self.tokens[self.pos]
        ):
            node = _Node("identifier", value=self.tokens[self.pos])
            self.pos += 1
        else:
            raise ValueError("Invalid KEGG module definition")
        return _Node("optional", [node]) if optional else node

    def _take(self, token: str) -> bool:
        if self.pos < len(self.tokens) and self.tokens[self.pos] == token:
            self.pos += 1
            return True
        return False


def ast_rows(module_id: str, root: _Node) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    counter = 0

    def visit(node: _Node, parent: int | None, position: int) -> None:
        nonlocal counter
        counter += 1
        node_id = counter
        rows.append(
            {
                "module_id": module_id,
                "node_id": node_id,
                "parent_node_id": parent,
                "position": position,
                "node_kind": node.kind,
                "member_namespace": None
                if node.value is None
                else ("ko" if node.value.startswith("K") else "module"),
                "member_id": node.value,
            }
        )
        for index, child in enumerate(node.children, 1):
            visit(child, node_id, index)

    visit(root, None, 1)
    return rows


def _empty(columns: Mapping[str, Any]) -> pl.DataFrame:
    return pl.DataFrame(schema=columns)


def _query_frame(connection: duckdb.DuckDBPyConnection, query: str) -> pl.DataFrame:
    cursor = connection.execute(query)
    columns = [item[0] for item in cursor.description]
    return pl.DataFrame(cursor.fetchall(), schema=columns, orient="row")


def write_duckdb(
    snapshot: MetabolicSnapshot,
    path: Path,
    *,
    if_exists: Literal["fail", "replace"] = "fail",
    include_source_hashes: bool = False,
) -> DuckDBWriteResult:
    from .publication import publish

    return publish(
        snapshot,
        path,
        if_exists=if_exists,
        include_source_hashes=include_source_hashes,
    )


def open_publication(path: Path) -> MetabolicPublication:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        con = duckdb.connect(str(path), read_only=True)
    except duckdb.Error as error:
        raise ValueError(f"Cannot open KEGG DuckDB publication: {path}") from error
    try:
        meta_tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM information_schema.tables WHERE table_schema='_bioextract'"
            ).fetchall()
        }
        required = {
            "metadata",
            "source_file",
            "table_info",
            "column_mapping",
            "validation_issue",
        }
        if not required <= meta_tables:
            missing = sorted(required - meta_tables)
            raise ValueError(
                f"DuckDB file is missing bioextract metadata tables: {missing}"
            )
        metadata: dict[str, str] = {
            str(key): str(value)
            for key, value in con.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
        }
        if metadata.get("bioextract.metadata_schema_version") != "1":
            raise ValueError("Unsupported KEGG metadata schema version")
        validate_duckdb_metadata_v1(con, metadata)
        if metadata.get("bioextract.source_schema_profile") != SOURCE_SCHEMA_PROFILE:
            raise ValueError("Unsupported KEGG source schema profile")
        if (
            metadata.get("bioextract.resource_name") != "kegg"
            or metadata.get("bioextract.scope") != "metabolic"
        ):
            raise ValueError(
                "DuckDB file is not a bioextract KEGG metabolic publication"
            )
        if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported KEGG metabolic resource schema version")
        tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM information_schema.tables WHERE table_schema='main'"
            ).fetchall()
        }
        recorded = dict(
            con.execute(
                "SELECT table_name, row_count FROM _bioextract.table_info"
            ).fetchall()
        )
        if tables != set(recorded):
            raise ValueError("KEGG DuckDB table inventory does not match metadata")
        for table, count in recorded.items():
            count_row = con.execute(f'SELECT count(*) FROM "{table}"').fetchone()
            if count_row is None or count_row[0] != count:
                raise ValueError(f"KEGG DuckDB row-count drift: {table}")
        capabilities: frozenset[str] = frozenset(
            value
            for value in metadata.get("bioextract.capabilities", "").split(",")
            if value
        )
        if metadata.get("bioextract.metadata_schema_version") == "1":
            if "bioextract.capabilities" not in metadata:
                raise ValueError(
                    "KEGG metadata v2 publication is missing capability metadata"
                )
            unknown = sorted(capabilities - set(_CAPABILITY_TABLES))
            if unknown:
                raise ValueError(
                    f"KEGG publication has unknown capabilities: {unknown}"
                )
            for capability in capabilities:
                if not (tables & _CAPABILITY_TABLES[capability]):
                    raise ValueError(
                        "KEGG capability metadata disagrees with table inventory: "
                        f"{capability}"
                    )
            for table in tables:
                if not any(
                    table in _CAPABILITY_TABLES[capability]
                    for capability in capabilities
                ):
                    raise ValueError(
                        "KEGG table inventory lacks a recorded capability owner: "
                        f"{table}"
                    )
        return MetabolicPublication(
            path,
            metadata,
            frozenset(tables),
            capabilities,
        )
    finally:
        con.close()


@dataclass(slots=True)
class KEGGMetabolicSelection:
    """Deferred reaction-centered selection over a metabolic publication.

    Instances are created by :meth:`KEGGDatabase.select_ids` or
    :meth:`KEGGDatabase.select_groups`.

    Examples:
        >>> selection = db.select_ids(  # doctest: +SKIP
        ...     ["CHEBI:15377"], namespace="chebi"
        ... )
        >>> selection.reactions().collect().height >= 1  # doctest: +SKIP
        True
    """

    publication: MetabolicPublication
    input_ids: tuple[str, ...]
    group_membership: tuple[tuple[str, str], ...] | None
    group_ids: tuple[str, ...]
    namespace: KEGGMetabolicNamespace
    include_obsolete: bool = False
    _matches_unique: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _matches: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _lineage: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _unmatched_unique: pl.DataFrame | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_ids(
        cls,
        publication: MetabolicPublication,
        ids: Iterable[str],
        namespace: KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KEGGMetabolicSelection:
        """Create a selection whose normalized identifiers are globally unique.

        Examples:
            Collapse prefixed and unprefixed forms before querying:

            >>> selection = KEGGMetabolicSelection.from_ids(  # doctest: +SKIP
            ...     publication,
            ...     ["cpd:C00001", "C00001"],
            ...     "kegg_compound",
            ... )
            >>> selection.input_ids  # doctest: +SKIP
            ('C00001',)
        """
        input_ids = tuple(
            sorted(
                {
                    normalized
                    for value in ids
                    if (normalized := _normalize_input(str(value), namespace))
                }
            )
        )
        return cls(
            publication=publication,
            input_ids=input_ids,
            group_membership=None,
            group_ids=(),
            namespace=namespace,
            include_obsolete=include_obsolete,
        )

    @classmethod
    def from_groups(
        cls,
        publication: MetabolicPublication,
        ids_by_group: Mapping[str, Iterable[str]],
        namespace: KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KEGGMetabolicSelection:
        """Create a grouped selection with one lookup row per normalized ID.

        Examples:
            Retain both group memberships while deduplicating the lookup ID:

            >>> selection = KEGGMetabolicSelection.from_groups(  # doctest: +SKIP
            ...     publication,
            ...     {"case": ["cpd:C00001"], "control": ["C00001"]},
            ...     "kegg_compound",
            ... )
            >>> selection.input_ids  # doctest: +SKIP
            ('C00001',)
        """
        group_ids = [str(group_id).strip() for group_id in ids_by_group]
        validate_group_ids(group_ids)
        membership = {
            (group_id, normalized)
            for group_id, values in zip(group_ids, ids_by_group.values(), strict=True)
            for value in values
            if (normalized := _normalize_input(str(value), namespace))
        }
        return cls(
            publication=publication,
            input_ids=tuple(sorted({input_id for _, input_id in membership})),
            group_membership=tuple(sorted(membership)),
            group_ids=tuple(sorted(group_ids)),
            namespace=namespace,
            include_obsolete=include_obsolete,
        )

    @property
    def is_grouped(self) -> bool:
        """Report whether extracted rows retain caller group labels.

        Examples:
            >>> selection.is_grouped  # doctest: +SKIP
            True
        """
        return self.group_membership is not None

    def _require(self, *tables: str) -> None:
        missing = [x for x in tables if x not in self.publication.tables]
        if missing:
            raise CapabilityError(
                f"KEGG metabolic publication lacks required relations: {missing}"
            )

    def matches(self) -> pl.LazyFrame:
        """Return canonical anchors resolved from every input identifier lazily.

        Examples:
            >>> selection.matches().select(  # doctest: +SKIP
            ...     "input_id", "entity_id"
            ... ).collect_schema().names()
            ['input_id', 'entity_id']
        """
        snapshot = copy.copy(self)
        prefix = ["group_id"] if self.is_grouped else []
        return register_deferred_frame_source(
            schema=_selection_schema(
                prefix
                + [
                    "input_id",
                    "input_namespace",
                    "match_type",
                    "entity_type",
                    "entity_id",
                ]
            ),
            frame=lambda: snapshot._eager_matches(),
        )

    def _eager_matches(self) -> pl.DataFrame:
        if self._matches is not None:
            return self._matches
        matches = self._resolve_unique_matches()
        self._matches = self._expand_groups(matches)
        return self._matches

    def _resolve_unique_matches(self) -> pl.DataFrame:
        if self._matches_unique is not None:
            return self._matches_unique
        ns = self.namespace
        if ns not in NAMESPACES:
            raise ValueError(
                f"Unknown KEGG metabolic namespace {ns!r}; available: {
                    list(NAMESPACES)
                }"
            )
        columns = pl.Schema(
            {
                "input_id": pl.String,
                "input_namespace": pl.String,
                "match_type": pl.String,
                "entity_type": pl.String,
                "entity_id": pl.String,
            }
        )
        if not self.input_ids:
            self._matches_unique = _empty(columns)
            return self._matches_unique

        with duckdb.connect(str(self.publication.path), read_only=True) as connection:
            _create_selection_input_table(connection, self.input_ids)
            self._matches_unique = (
                self._query_unique_matches(connection)
                .cast(columns)
                .unique()
                .sort(list(columns))
            )
        return self._matches_unique

    def _query_unique_matches(
        self, connection: duckdb.DuckDBPyConnection
    ) -> pl.DataFrame:
        ns = self.namespace
        direct = {
            "kegg_compound": ("compound", "compound_id", "compound"),
            "kegg_reaction": ("reaction", "reaction_id", "reaction"),
            "ec": ("enzyme", "ec_number", "enzyme"),
            "kegg_module": ("module", "module_id", "module"),
        }
        xref = {
            "chebi": (
                "compound_cross_reference",
                "external_id",
                "compound_id",
                "compound",
            ),
            "rhea": (
                "reaction_cross_reference",
                "external_id",
                "reaction_id",
                "reaction",
            ),
            "pubchem": ("compound_pubchem", "pubchem_id", "compound_id", "compound"),
        }
        relation = {
            "ko": ("reaction_ko", "ko_id", "reaction_id", "reaction"),
            "kegg_pathway": (
                "reaction_pathway",
                "pathway_id",
                "reaction_id",
                "reaction",
            ),
        }
        table: str
        column: str
        target: str
        entity: str
        namespace_filter = ""
        if ns == "kegg_compound" and "compound" not in self.publication.tables:
            if "compound_reaction" in self.publication.tables:
                table, column, target, entity = (
                    "compound_reaction",
                    "compound_id",
                    "compound_id",
                    "compound",
                )
            else:
                self._require("reaction_participant")
                table, column, target, entity = (
                    "reaction_participant",
                    "participant_id",
                    "participant_id",
                    "compound",
                )
                namespace_filter = (
                    " AND target_row.participant_namespace='kegg_compound'"
                )
        elif ns == "ec" and "enzyme" not in self.publication.tables:
            self._require("reaction_enzyme")
            table, column, target, entity = (
                "reaction_enzyme",
                "ec_number",
                "reaction_id",
                "reaction",
            )
        elif ns == "kegg_module" and "module" not in self.publication.tables:
            self._require("reaction_module")
            table, column, target, entity = (
                "reaction_module",
                "module_id",
                "reaction_id",
                "reaction",
            )
        elif ns == "ec":
            self._require("enzyme")
            return self._query_unique_ec_matches(connection)
        elif ns in direct:
            table, column, entity = direct[ns]
            target = column
            self._require(table)
        elif ns in xref:
            table, column, target, entity = xref[ns]
            self._require(table)
            if ns in {"chebi", "rhea"}:
                namespace_filter = f" AND target_row.namespace='{ns}'"
        else:
            table, column, target, entity = relation[ns]
            self._require(table)

        return _query_frame(
            connection,
            f"""
            SELECT DISTINCT
                input.input_id AS "input_id",
                '{ns}' AS "input_namespace",
                'exact' AS "match_type",
                '{entity}' AS entity_type,
                target_row.{target} AS entity_id
            FROM _input_id AS input
            JOIN {table} AS target_row
              ON target_row.{column} = input.input_id
            WHERE true{namespace_filter}
            ORDER BY input_id, entity_id
            """,
        )

    def _query_unique_ec_matches(
        self, connection: duckdb.DuckDBPyConnection
    ) -> pl.DataFrame:
        if self.include_obsolete:
            return _query_frame(
                connection,
                """
                SELECT DISTINCT
                    input.input_id AS "input_id",
                    'ec' AS "input_namespace",
                    'exact' AS "match_type",
                    'enzyme' AS entity_type,
                    enzyme.ec_number AS entity_id
                FROM _input_id AS input
                JOIN enzyme ON enzyme.ec_number = input.input_id
                ORDER BY input_id, entity_id
                """,
            )

        replacement = (
            """
            UNION ALL
            SELECT DISTINCT
                path.input_id AS "input_id",
                'ec' AS "input_namespace",
                'replacement' AS "match_type",
                'enzyme' AS entity_type,
                enzyme.ec_number AS entity_id
            FROM replacement_path AS path
            JOIN enzyme ON enzyme.ec_number = path.ec_number
            WHERE enzyme.status = 'active'
            """
            if "enzyme_replacement" in self.publication.tables
            else ""
        )
        recursive_cte = (
            """
            WITH RECURSIVE replacement_path(input_id, ec_number) AS (
                SELECT input.input_id, edge.replacement_ec_number
                FROM _input_id AS input
                JOIN enzyme_replacement AS edge
                  ON edge.ec_number = input.input_id
                UNION
                SELECT path.input_id, edge.replacement_ec_number
                FROM replacement_path AS path
                JOIN enzyme_replacement AS edge
                  ON edge.ec_number = path.ec_number
            )
            """
            if replacement
            else ""
        )
        return _query_frame(
            connection,
            f"""
            {recursive_cte}
            SELECT DISTINCT
                input.input_id AS "input_id",
                'ec' AS "input_namespace",
                'exact' AS "match_type",
                'enzyme' AS entity_type,
                enzyme.ec_number AS entity_id
            FROM _input_id AS input
            JOIN enzyme ON enzyme.ec_number = input.input_id
            WHERE enzyme.status = 'active'
            {replacement}
            ORDER BY input_id, match_type, entity_id
            """,
        )

    def _expand_groups(self, frame: pl.DataFrame) -> pl.DataFrame:
        if self.group_membership is None:
            return frame
        membership = pl.DataFrame(
            self.group_membership,
            schema={"group_id": pl.String, "input_id": pl.String},
            orient="row",
        )
        return (
            membership.join(frame, on="input_id", how="inner")
            .select("group_id", *frame.columns)
            .sort("group_id", *frame.columns)
        )

    def _reaction_lineage(self) -> pl.DataFrame:
        if self._lineage is not None:
            return self._lineage
        matches = self._eager_matches()
        prefix = ["group_id"] if self.is_grouped else []
        schema = pl.Schema(
            dict.fromkeys(
                (
                    *prefix,
                    "input_id",
                    "input_namespace",
                    "anchor_type",
                    "anchor_id",
                    "reaction_id",
                ),
                pl.String,
            )
        )
        if matches.is_empty():
            self._lineage = _empty(schema)
            return self._lineage

        entity_types = set(matches["entity_type"].to_list())
        clauses = [
            """
            SELECT
                selected.group_id,
                selected.input_id,
                selected.input_namespace,
                selected.entity_type AS anchor_type,
                selected.entity_id AS anchor_id,
                selected.entity_id AS reaction_id
            FROM _selected_anchor AS selected
            WHERE selected.entity_type = 'reaction'
            """
        ]
        if "compound" in entity_types:
            if "compound_reaction" in self.publication.tables:
                clauses.append(
                    """
                    SELECT
                        selected.group_id,
                        selected.input_id,
                        selected.input_namespace,
                        selected.entity_type,
                        selected.entity_id,
                        relation.reaction_id
                    FROM _selected_anchor AS selected
                    JOIN compound_reaction AS relation
                      ON relation.compound_id = selected.entity_id
                    WHERE selected.entity_type = 'compound'
                    """
                )
            else:
                self._require("reaction_participant")
                clauses.append(
                    """
                    SELECT
                        selected.group_id,
                        selected.input_id,
                        selected.input_namespace,
                        selected.entity_type,
                        selected.entity_id,
                        participant.reaction_id
                    FROM _selected_anchor AS selected
                    JOIN reaction_participant AS participant
                      ON participant.participant_id = selected.entity_id
                     AND participant.participant_namespace = 'kegg_compound'
                    WHERE selected.entity_type = 'compound'
                    """
                )
        if "enzyme" in entity_types:
            self._require("reaction_enzyme")
            clauses.append(
                """
                SELECT
                    selected.group_id,
                    selected.input_id,
                    selected.input_namespace,
                    selected.entity_type,
                    selected.entity_id,
                    relation.reaction_id
                FROM _selected_anchor AS selected
                JOIN reaction_enzyme AS relation
                  ON relation.ec_number = selected.entity_id
                WHERE selected.entity_type = 'enzyme'
                """
            )
        if "module" in entity_types:
            self._require("reaction_module")
            clauses.append(
                """
                SELECT
                    selected.group_id,
                    selected.input_id,
                    selected.input_namespace,
                    selected.entity_type,
                    selected.entity_id,
                    relation.reaction_id
                FROM _selected_anchor AS selected
                JOIN reaction_module AS relation
                  ON relation.module_id = selected.entity_id
                WHERE selected.entity_type = 'module'
                """
            )

        with duckdb.connect(str(self.publication.path), read_only=True) as connection:
            _create_selected_anchor_table(connection, matches)
            lineage = _query_frame(
                connection,
                f"""
                SELECT DISTINCT
                    group_id AS "group_id",
                    input_id AS "input_id",
                    input_namespace AS "input_namespace",
                    anchor_type AS anchor_type,
                    anchor_id AS anchor_id,
                    reaction_id AS reaction_id
                FROM ({" UNION ALL ".join(clauses)})
                ORDER BY group_id, input_id, anchor_type, anchor_id, reaction_id
                """,
            )
        self._lineage = (
            lineage.select(*schema).cast(schema)
            if self.is_grouped
            else lineage.drop("group_id").select(*schema).cast(schema)
        )
        return self._lineage

    def _extract(self, table: str, join_column: str) -> pl.DataFrame:
        self._require(table)
        lineage = self._reaction_lineage()
        if lineage.is_empty():
            return lineage
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            frame = _query_frame(con, f"SELECT * FROM {table}")
        result = lineage.join(
            frame, left_on="reaction_id", right_on=join_column, how="inner"
        )
        return result

    def reactions(self) -> pl.LazyFrame:
        """Return selected reactions with input and anchor lineage lazily.

        Examples:
            >>> selection.reactions().select(  # doctest: +SKIP
            ...     "reaction_id", "equation"
            ... ).collect_schema().names()
            ['reaction_id', 'equation']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction")
            ),
            frame=lambda: snapshot._eager_reactions(),
        )

    def _eager_reactions(self) -> pl.DataFrame:
        return self._extract("reaction", "reaction_id")

    def participants(self) -> pl.LazyFrame:
        """Return ordered left/right participants lazily.

        Examples:
            >>> selection.participants().select(  # doctest: +SKIP
            ...     "reaction_id", "side", "participant_id"
            ... ).collect_schema().names()
            ['reaction_id', 'side', 'participant_id']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction_participant")
            ),
            frame=lambda: snapshot._eager_participants(),
        )

    def _eager_participants(self) -> pl.DataFrame:
        return self._extract("reaction_participant", "reaction_id")

    def enzymes(self) -> pl.LazyFrame:
        """Return EC links owned by the selected reactions lazily.

        Examples:
            >>> selection.enzymes().select(  # doctest: +SKIP
            ...     "reaction_id", "ec_number"
            ... ).collect_schema().names()
            ['reaction_id', 'ec_number']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction_enzyme")
            ),
            frame=lambda: snapshot._eager_enzymes(),
        )

    def _eager_enzymes(self) -> pl.DataFrame:
        return self._extract("reaction_enzyme", "reaction_id")

    def kos(self) -> pl.LazyFrame:
        """Return KO links owned by the selected reactions lazily.

        Examples:
            >>> selection.kos().select(  # doctest: +SKIP
            ...     "reaction_id", "ko_id"
            ... ).collect_schema().names()
            ['reaction_id', 'ko_id']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction_ko")
            ),
            frame=lambda: snapshot._eager_kos(),
        )

    def _eager_kos(self) -> pl.DataFrame:
        return self._extract("reaction_ko", "reaction_id")

    def modules(self) -> pl.LazyFrame:
        """Return module memberships of the selected reactions lazily.

        Examples:
            >>> selection.modules().select(  # doctest: +SKIP
            ...     "reaction_id", "module_id"
            ... ).collect_schema().names()
            ['reaction_id', 'module_id']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction_module")
            ),
            frame=lambda: snapshot._eager_modules(),
        )

    def _eager_modules(self) -> pl.DataFrame:
        return self._extract("reaction_module", "reaction_id")

    def compounds(self) -> pl.LazyFrame:
        """Return compound participants and facts lazily.

        Examples:
            >>> selection.compounds().select(  # doctest: +SKIP
            ...     "participant_id", "name"
            ... ).collect_schema().names()
            ['participant_id', 'name']
        """
        snapshot = copy.copy(self)
        prefix = ["group_id"] if self.is_grouped else []
        participant_columns = (
            prefix
            + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
            + _table_columns(self.publication, "reaction_participant")
        )
        columns = participant_columns + [
            column
            for column in _table_columns(self.publication, "compound")
            if column not in {"compound_id"}
        ]
        return register_deferred_frame_source(
            schema=_selection_schema(columns),
            frame=lambda: snapshot._eager_compounds(),
        )

    def _eager_compounds(self) -> pl.DataFrame:
        participants = self._eager_participants().filter(
            pl.col("participant_namespace") == "kegg_compound"
        )
        self._require("compound")
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            compounds = _query_frame(con, "SELECT * FROM compound")
        return participants.join(
            compounds,
            left_on="participant_id",
            right_on="compound_id",
            how="left",
        )

    def pathway_memberships(self) -> pl.LazyFrame:
        """Return reference-pathway memberships lazily.

        Examples:
            >>> selection.pathway_memberships().select(  # doctest: +SKIP
            ...     "reaction_id", "pathway_id"
            ... ).collect_schema().names()
            ['reaction_id', 'pathway_id']
        """
        snapshot = copy.copy(self)
        return register_deferred_frame_source(
            schema=_selection_schema(
                (["group_id"] if self.is_grouped else [])
                + ["input_id", "input_namespace", "anchor_type", "anchor_id"]
                + _table_columns(self.publication, "reaction_pathway")
            ),
            frame=lambda: snapshot._eager_pathway_memberships(),
        )

    def _eager_pathway_memberships(self) -> pl.DataFrame:
        return self._extract("reaction_pathway", "reaction_id")

    def cross_references(self) -> pl.LazyFrame:
        """Return compound and reaction cross-references lazily.

        Examples:
            >>> selection.cross_references().select(  # doctest: +SKIP
            ...     "namespace", "external_id"
            ... ).collect_schema().names()
            ['namespace', 'external_id']
        """
        snapshot = copy.copy(self)
        columns = ["group_id"] if self.is_grouped else []
        columns += ["input_id", "input_namespace", "anchor_type", "anchor_id"]
        columns += _table_columns(self.publication, "reaction_participant")
        if "compound" in self.publication.tables:
            columns += [
                column
                for column in _table_columns(self.publication, "compound")
                if column != "compound_id"
            ]
        for table in ("compound_cross_reference", "reaction_cross_reference"):
            if table in self.publication.tables:
                columns += [
                    column
                    for column in _table_columns(self.publication, table)
                    if column not in columns
                ]
        return register_deferred_frame_source(
            schema=_selection_schema(dict.fromkeys(columns)),
            frame=lambda: snapshot._eager_cross_references(),
        )

    def _eager_cross_references(self) -> pl.DataFrame:
        frames: list[pl.DataFrame] = []
        compounds = self._eager_compounds()
        reactions = self._eager_reactions()
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            if "compound_cross_reference" in self.publication.tables:
                frames.append(
                    compounds.join(
                        _query_frame(con, "SELECT * FROM compound_cross_reference"),
                        left_on="participant_id",
                        right_on="compound_id",
                        how="inner",
                    )
                )
            if "reaction_cross_reference" in self.publication.tables:
                frames.append(
                    reactions.join(
                        _query_frame(con, "SELECT * FROM reaction_cross_reference"),
                        left_on="reaction_id",
                        right_on="reaction_id",
                        how="inner",
                    )
                )
        return pl.concat(frames, how="diagonal_relaxed") if frames else pl.DataFrame()

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return normalized inputs that resolved to no canonical anchor lazily.

        Examples:
            >>> selection.unmatched_ids().select(  # doctest: +SKIP
            ...     "input_id", "reason"
            ... ).collect_schema().names()
            ['input_id', 'reason']
        """
        snapshot = copy.copy(self)
        prefix = ["group_id"] if self.is_grouped else []
        return register_deferred_frame_source(
            schema=_selection_schema(prefix + ["input_id", "reason"]),
            frame=lambda: snapshot._eager_unmatched_ids(),
        )

    def _eager_unmatched_ids(self) -> pl.DataFrame:
        return self._expand_groups(self._resolve_unique_unmatched())

    def _resolve_unique_unmatched(self) -> pl.DataFrame:
        if self._unmatched_unique is not None:
            return self._unmatched_unique
        matched = set(self._resolve_unique_matches()["input_id"].to_list())
        missing = tuple(
            input_id for input_id in self.input_ids if input_id not in matched
        )
        schema = {"input_id": pl.String, "reason": pl.String}
        if not missing:
            self._unmatched_unique = _empty(schema)
            return self._unmatched_unique

        reasons = dict.fromkeys(missing, "not_found")
        if self.namespace == "ec" and "enzyme" in self.publication.tables:
            with duckdb.connect(
                str(self.publication.path), read_only=True
            ) as connection:
                _create_selection_input_table(connection, missing)
                statuses = connection.execute(
                    """
                    SELECT input.input_id, enzyme.status
                    FROM _input_id AS input
                    JOIN enzyme ON enzyme.ec_number = input.input_id
                    """
                ).fetchall()
            for input_id, status in statuses:
                if status == "deleted" and not self.include_obsolete:
                    reasons[str(input_id)] = "obsolete_excluded"
                elif status == "transferred":
                    reasons[str(input_id)] = "invalid_canonical_target"

        self._unmatched_unique = pl.DataFrame(
            {
                "input_id": missing,
                "reason": tuple(reasons[input_id] for input_id in missing),
            },
            schema=schema,
        )
        return self._unmatched_unique


def _create_selection_input_table(
    connection: duckdb.DuckDBPyConnection, input_ids: Sequence[str]
) -> None:
    connection.execute("CREATE TEMP TABLE _input_id(input_id VARCHAR PRIMARY KEY)")
    connection.executemany(
        "INSERT INTO _input_id VALUES (?)",
        [(input_id,) for input_id in input_ids],
    )


def _create_selected_anchor_table(
    connection: duckdb.DuckDBPyConnection, matches: pl.DataFrame
) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _selected_anchor(
            group_id VARCHAR,
            input_id VARCHAR NOT NULL,
            input_namespace VARCHAR NOT NULL,
            entity_type VARCHAR NOT NULL,
            entity_id VARCHAR NOT NULL
        )
        """
    )
    connection.executemany(
        "INSERT INTO _selected_anchor VALUES (?, ?, ?, ?, ?)",
        [
            (
                row.get("group_id"),
                row["input_id"],
                row["input_namespace"],
                row["entity_type"],
                row["entity_id"],
            )
            for row in matches.iter_rows(named=True)
        ],
    )


def _normalize_input(value: str, namespace: str) -> str:
    value = str(value).strip()
    if not value:
        return ""
    prefixes = {
        "kegg_compound": "cpd:",
        "kegg_reaction": "rn:",
        "kegg_module": "md:",
        "kegg_pathway": "path:",
    }
    value = value.removeprefix(prefixes.get(namespace, ""))
    if not value:
        return ""
    if namespace == "chebi":
        identifier = value.upper().removeprefix("CHEBI:")
        value = f"CHEBI:{identifier}" if identifier else ""
    if namespace == "rhea":
        identifier = value.upper().removeprefix("RHEA:")
        value = f"RHEA:{identifier}" if identifier else ""
    return value


def validate_selection_namespace(
    publication: MetabolicPublication,
    namespace: KEGGMetabolicNamespace,
) -> None:
    tables = publication.tables
    cross_reference_namespaces: set[str] = set()
    with duckdb.connect(str(publication.path), read_only=True) as connection:
        for table in ("compound_cross_reference", "reaction_cross_reference"):
            if table in tables:
                cross_reference_namespaces.update(
                    str(row[0])
                    for row in connection.execute(
                        f"SELECT DISTINCT namespace FROM {table}"
                    ).fetchall()
                )
    available = {
        *(
            ("kegg_compound",)
            if tables & {"compound", "compound_reaction", "reaction_participant"}
            else ()
        ),
        *(("chebi",) if "chebi" in cross_reference_namespaces else ()),
        *(("pubchem",) if "compound_pubchem" in tables else ()),
        *(("kegg_reaction",) if "reaction" in tables else ()),
        *(("rhea",) if "rhea" in cross_reference_namespaces else ()),
        *(("ec",) if tables & {"enzyme", "reaction_enzyme"} else ()),
        *(("ko",) if "reaction_ko" in tables else ()),
        *(("kegg_module",) if tables & {"module", "reaction_module"} else ()),
        *(("kegg_pathway",) if "reaction_pathway" in tables else ()),
    }
    if namespace not in available:
        raise CapabilityError(
            f"KEGG metabolic namespace {namespace!r} is unavailable; "
            f"available namespaces: {sorted(available)}"
        )


def evaluate_modules(
    publication: MetabolicPublication, ko_ids: Iterable[str]
) -> pl.DataFrame:
    if "module_definition_node" not in publication.tables:
        raise CapabilityError("KEGG metabolic publication lacks module definitions")
    available = {str(x).strip().removeprefix("ko:") for x in ko_ids}
    with duckdb.connect(str(publication.path), read_only=True) as con:
        rows = con.execute(
            "SELECT module_id,node_id,parent_node_id,position,node_kind,member_namespace,member_id FROM module_definition_node ORDER BY module_id,node_id"
        ).fetchall()
    by_module: dict[str, dict[int, tuple[int | None, str, str | None]]] = {}
    for module, node, parent, _position, kind, _ns, member in rows:
        by_module.setdefault(module, {})[node] = (parent, kind, member)
    roots: dict[str, int] = {}
    children_by_module: dict[str, dict[int, list[int]]] = {}
    for module, nodes in by_module.items():
        children: dict[int, list[int]] = {}
        for node, (parent, _, _) in nodes.items():
            if parent is None:
                roots[module] = node
            else:
                children.setdefault(parent, []).append(node)
        children_by_module[module] = children

    def identifier_satisfied(member: str | None, stack: tuple[str, ...]) -> bool:
        if member is None:
            return False
        if member.startswith("K"):
            return member in available
        return module_satisfied(member, stack)

    def module_satisfied(module_id: str, stack: tuple[str, ...]) -> bool:
        if module_id in stack:
            cycle = " -> ".join((*stack, module_id))
            raise ValueError(f"Cyclic KEGG module definition reference: {cycle}")
        nodes = by_module.get(module_id)
        root = roots.get(module_id)
        if nodes is None or root is None:
            return False
        children = children_by_module[module_id]

        def satisfied(node: int) -> bool:
            _, kind, member = nodes[node]
            values = [satisfied(child) for child in children.get(node, [])]
            if kind == "identifier":
                return identifier_satisfied(member, (*stack, module_id))
            if kind == "alternative":
                return any(values)
            if kind == "optional":
                return True
            return all(values)

        return satisfied(root)

    result: list[dict[str, Any]] = []
    for module, nodes in by_module.items():
        children = children_by_module[module]
        root = roots[module]

        def block_satisfied(
            node: int,
            module_nodes: Mapping[int, tuple[int | None, str, str | None]] = nodes,
            node_children: Mapping[int, list[int]] = children,
            module_id: str = module,
        ) -> bool:
            _, kind, member = module_nodes[node]
            values = [block_satisfied(child) for child in node_children.get(node, [])]
            if kind == "identifier":
                return identifier_satisfied(member, (module_id,))
            if kind == "alternative":
                return any(values)
            if kind == "optional":
                return True
            return all(values)

        blocks = (
            [
                node
                for node in children.get(root, [root])
                if nodes[node][1] != "optional"
            ]
            if nodes[root][1] == "sequence"
            else ([] if nodes[root][1] == "optional" else [root])
        )
        hits = [block_satisfied(node) for node in blocks]
        result.append(
            {
                "module_id": module,
                "required_block_count": len(blocks),
                "satisfied_block_count": sum(hits),
                "is_complete": all(hits),
                "missing_block_indexes": [
                    i for i, hit in enumerate(hits, 1) if not hit
                ],
            }
        )
    return pl.DataFrame(result).sort("module_id")
