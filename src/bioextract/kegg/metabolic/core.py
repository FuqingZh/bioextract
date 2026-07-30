from __future__ import annotations

import gzip
import os
import re
import tarfile
import zipfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import duckdb
import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    validate_duckdb_metadata_v3,
    validate_duckdb_validation_state,
)

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


class KEGGMetabolicCapabilityError(RuntimeError):
    """Raised when a partial KEGG publication lacks a required relation.

    Examples:
        >>> error = KEGGMetabolicCapabilityError("reaction_ko is unavailable")
        >>> "reaction_ko" in str(error)
        True
    """


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
    entries: bool = False,
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
            pattern = "*.keg" if entries else "*"
            result.extend(sorted(p for p in path.rglob(pattern) if p.is_file()))
        else:
            result.append(path)
    return tuple(sorted(result))


def from_metabolic_files(
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
    values: dict[str, Any] = locals()
    sources = {
        role: _paths(values[role], entries=role.endswith("_entries"))
        for role in (
            *[f"{x}_list" for x in ENTRY_ROLES],
            *[f"{x}_entries" for x in ENTRY_ROLES],
            *RELATION_ROLES,
        )
        if values[role] is not None
    }
    if not sources:
        raise ValueError("At least one KEGG metabolic input must be provided")
    return MetabolicSnapshot(sources=sources, release_version=release_version)


def from_metabolic_release(
    source: os.PathLike[str] | str,
    *,
    release_version: str | None = None,
) -> MetabolicSnapshot:
    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_file():
        if not (zipfile.is_zipfile(path) or tarfile.is_tarfile(path)):
            raise ValueError(f"Unsupported KEGG metabolic release archive: {path}")
        return MetabolicSnapshot(
            sources={},
            release_version=release_version,
            complete_release=True,
            archive=path,
        )
    root = path / "raw" if (path / "raw").is_dir() else path
    sources = discover_release_layout(root)
    missing = [
        role
        for role in (*[f"{x}_entries" for x in ENTRY_ROLES], *RELATION_ROLES)
        if role not in sources
    ]
    if missing:
        raise ValueError(f"Incomplete KEGG metabolic release; missing roles: {missing}")
    return MetabolicSnapshot(
        sources=sources,
        release_version=release_version,
        complete_release=True,
    )


def discover_release_layout(root: Path) -> dict[str, tuple[Path, ...]]:
    candidates = [root]
    if (root / "raw").is_dir():
        candidates.insert(0, root / "raw")
    candidates.extend(sorted(path for path in root.rglob("raw") if path.is_dir()))
    raw = next(
        (
            candidate
            for candidate in candidates
            if all((candidate / family / "entries").is_dir() for family in ENTRY_ROLES)
        ),
        root,
    )
    found: dict[str, tuple[Path, ...]] = {}
    for family in ENTRY_ROLES:
        file_list = raw / family / "list.tsv"
        dir_entries = raw / family / "entries"
        if file_list.is_file():
            found[f"{family}_list"] = (file_list,)
        if dir_entries.is_dir():
            entries = tuple(sorted(dir_entries.glob("*.keg")))
            if entries:
                found[f"{family}_entries"] = entries
    for role in RELATION_ROLES:
        relation = raw / f"{role}.tsv"
        if relation.is_file():
            found[role] = (relation,)
    return found


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


def _public_frame(frame: pl.DataFrame) -> pl.DataFrame:
    renamed = {
        column: "".join(part.capitalize() for part in column.split("_"))
        for column in frame.columns
        if "_" in column
    }
    renamed.update({"ec_number": "EcNumber", "ko_id": "KoId"})
    return frame.rename({key: value for key, value in renamed.items() if key in frame})


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
        required = {"metadata", "source_file", "table_info", "column_mapping"}
        if not required <= meta_tables:
            raise ValueError("DuckDB file is missing bioextract metadata tables")
        metadata: dict[str, str] = {
            str(key): str(value)
            for key, value in con.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
        }
        if metadata.get("bioextract.metadata_schema_version") not in {"1", "2", "3"}:
            raise ValueError("Unsupported KEGG metadata schema version")
        if metadata.get("bioextract.metadata_schema_version") == "3":
            validate_duckdb_metadata_v3(con, metadata)
            required_v3 = {
                "bioextract.resource_name",
                "bioextract.resource_schema_version",
                "bioextract.source_schema_profile",
                "bioextract.package_version",
                "bioextract.generated_at",
                "bioextract.validation_status",
                "bioextract.validation_issue_count",
                "bioextract.sources",
            }
            missing_v3 = sorted(required_v3 - set(metadata))
            if missing_v3:
                raise ValueError(f"KEGG metadata v3 is missing keys: {missing_v3}")
            if (
                metadata.get("bioextract.source_schema_profile")
                != SOURCE_SCHEMA_PROFILE
            ):
                raise ValueError("Unsupported KEGG source schema profile")
        if metadata.get("bioextract.metadata_schema_version") == "2":
            validate_duckdb_validation_state(con, metadata)
        if (
            metadata.get("bioextract.resource_name") != "kegg"
            or metadata.get("bioextract.scope") != "metabolic"
        ):
            raise ValueError(
                "DuckDB file is not a bioextract KEGG metabolic publication"
            )
        metadata_version = metadata.get("bioextract.metadata_schema_version")
        resource_schema_key = (
            "bioextract.resource_schema_version"
            if metadata_version == "3"
            else "bioextract.schema_version"
        )
        if metadata.get(resource_schema_key) != SCHEMA_VERSION:
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
        if metadata.get("bioextract.metadata_schema_version") in {"2", "3"}:
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
        >>> selection.extract_reactions().height >= 1  # doctest: +SKIP
        True
    """

    publication: MetabolicPublication
    inputs: tuple[tuple[str | None, str], ...]
    namespace: KEGGMetabolicNamespace
    include_obsolete: bool = False
    _matches: pl.DataFrame | None = field(default=None, init=False, repr=False)

    @property
    def is_grouped(self) -> bool:
        return any(group is not None for group, _ in self.inputs)

    def _require(self, *tables: str) -> None:
        missing = [x for x in tables if x not in self.publication.tables]
        if missing:
            raise KEGGMetabolicCapabilityError(
                f"KEGG metabolic publication lacks required relations: {missing}"
            )

    def extract_matches(self) -> pl.DataFrame:
        """Return canonical anchors resolved from every input identifier.

        Examples:
            >>> selection.extract_matches().select(  # doctest: +SKIP
            ...     "InputId", "EntityId"
            ... ).head(1)
        """
        if self._matches is not None:
            return self._matches
        ns = self.namespace
        if ns not in NAMESPACES:
            raise ValueError(
                f"Unknown KEGG metabolic namespace {ns!r}; available: {
                    list(NAMESPACES)
                }"
            )
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
        rows: list[dict[str, Any]] = []
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            for group, raw in self.inputs:
                value = _normalize_input(raw, ns)
                if ns == "kegg_compound" and "compound" not in self.publication.tables:
                    if "compound_reaction" in self.publication.tables:
                        table, column, target, entity = (
                            "compound_reaction",
                            "compound_id",
                            "compound_id",
                            "compound",
                        )
                    else:
                        table, column, target, entity = (
                            "reaction_participant",
                            "participant_id",
                            "participant_id",
                            "compound",
                        )
                    ids = [
                        row[0]
                        for row in con.execute(
                            f"SELECT DISTINCT {target} FROM {table} WHERE {column}=?",
                            [value],
                        ).fetchall()
                    ]
                elif ns == "ec" and "enzyme" not in self.publication.tables:
                    table, column, target, entity = (
                        "reaction_enzyme",
                        "ec_number",
                        "reaction_id",
                        "reaction",
                    )
                    ids = [
                        row[0]
                        for row in con.execute(
                            f"SELECT DISTINCT {target} FROM {table} WHERE {column}=?",
                            [value],
                        ).fetchall()
                    ]
                elif ns == "kegg_module" and "module" not in self.publication.tables:
                    table, column, target, entity = (
                        "reaction_module",
                        "module_id",
                        "reaction_id",
                        "reaction",
                    )
                    ids = [
                        row[0]
                        for row in con.execute(
                            f"SELECT DISTINCT {target} FROM {table} WHERE {column}=?",
                            [value],
                        ).fetchall()
                    ]
                elif ns in direct:
                    table, column, entity = direct[ns]
                    self._require(table)
                    if ns == "ec":
                        status_clause = (
                            "" if self.include_obsolete else " AND status='active'"
                        )
                        ids = [
                            r[0]
                            for r in con.execute(
                                f"SELECT {column} FROM {table} "
                                f"WHERE {column}=?{status_clause}",
                                [value],
                            ).fetchall()
                        ]
                        if (
                            not ids
                            and not self.include_obsolete
                            and "enzyme_replacement" in self.publication.tables
                        ):
                            ids = [
                                r[0]
                                for r in con.execute(
                                    """
                                    WITH RECURSIVE replacement_path(ec_number) AS (
                                        SELECT replacement_ec_number
                                        FROM enzyme_replacement
                                        WHERE ec_number = ?
                                        UNION
                                        SELECT edge.replacement_ec_number
                                        FROM enzyme_replacement AS edge
                                        JOIN replacement_path AS path
                                          ON edge.ec_number = path.ec_number
                                    )
                                    SELECT DISTINCT enzyme.ec_number
                                    FROM replacement_path
                                    JOIN enzyme USING (ec_number)
                                    WHERE enzyme.status = 'active'
                                    ORDER BY enzyme.ec_number
                                    """,
                                    [value],
                                ).fetchall()
                            ]
                    else:
                        ids = [
                            r[0]
                            for r in con.execute(
                                f"SELECT {column} FROM {table} WHERE {column}=?",
                                [value],
                            ).fetchall()
                        ]
                elif ns in xref:
                    table, column, target, entity = xref[ns]
                    self._require(table)
                    ids = [
                        r[0]
                        for r in con.execute(
                            f"SELECT {target} FROM {table} WHERE {column}=?", [value]
                        ).fetchall()
                    ]
                else:
                    table, column, target, entity = relation[ns]
                    self._require(table)
                    ids = [
                        r[0]
                        for r in con.execute(
                            f"SELECT DISTINCT {target} FROM {table} WHERE {column}=?",
                            [value],
                        ).fetchall()
                    ]
                for ident in ids:
                    rows.append(
                        {
                            "GroupId": group,
                            "InputId": value,
                            "InputNamespace": ns,
                            "MatchType": (
                                "replacement"
                                if ns == "ec" and entity == "enzyme" and ident != value
                                else "exact"
                            ),
                            "EntityType": entity,
                            "EntityId": ident,
                        }
                    )
        columns = [
            "GroupId",
            "InputId",
            "InputNamespace",
            "MatchType",
            "EntityType",
            "EntityId",
        ]
        self._matches = (
            pl.DataFrame(rows, schema=dict.fromkeys(columns, pl.String))
            if rows
            else _empty(dict.fromkeys(columns, pl.String))
        )
        if not self.is_grouped:
            self._matches = self._matches.drop("GroupId")
        return self._matches

    def _reaction_lineage(self) -> pl.DataFrame:
        matches = self.extract_matches()
        prefix = ["GroupId"] if self.is_grouped else []
        rows: list[dict[str, Any]] = []
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            for match in matches.iter_rows(named=True):
                entity, ident = match["EntityType"], match["EntityId"]
                reactions: list[str]
                if entity == "reaction":
                    reactions = [ident]
                elif entity == "compound":
                    if "compound_reaction" in self.publication.tables:
                        query = (
                            "SELECT reaction_id FROM compound_reaction "
                            "WHERE compound_id=?"
                        )
                    else:
                        self._require("reaction_participant")
                        query = (
                            "SELECT reaction_id FROM reaction_participant "
                            "WHERE participant_namespace='kegg_compound' "
                            "AND participant_id=?"
                        )
                    reactions = [r[0] for r in con.execute(query, [ident]).fetchall()]
                elif entity == "enzyme":
                    self._require("reaction_enzyme")
                    reactions = [
                        r[0]
                        for r in con.execute(
                            "SELECT reaction_id FROM reaction_enzyme WHERE ec_number=?",
                            [ident],
                        ).fetchall()
                    ]
                elif entity == "module":
                    self._require("reaction_module")
                    reactions = [
                        r[0]
                        for r in con.execute(
                            "SELECT reaction_id FROM reaction_module WHERE module_id=?",
                            [ident],
                        ).fetchall()
                    ]
                else:
                    reactions = []
                for reaction_id in reactions:
                    rows.append(
                        {
                            **{
                                k: match[k]
                                for k in (*prefix, "InputId", "InputNamespace")
                            },
                            "AnchorType": entity,
                            "AnchorId": ident,
                            "ReactionId": reaction_id,
                        }
                    )
        schema = dict.fromkeys(
            (
                *prefix,
                "InputId",
                "InputNamespace",
                "AnchorType",
                "AnchorId",
                "ReactionId",
            ),
            pl.String,
        )
        return pl.DataFrame(rows, schema=schema).unique() if rows else _empty(schema)

    def _extract(
        self, table: str, join_column: str, rename: Mapping[str, str] | None = None
    ) -> pl.DataFrame:
        self._require(table)
        lineage = self._reaction_lineage()
        if lineage.is_empty():
            return lineage
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            frame = _query_frame(con, f"SELECT * FROM {table}")
        result = lineage.join(
            frame, left_on="ReactionId", right_on=join_column, how="inner"
        )
        return _public_frame(result.rename(rename or {}))

    def extract_reactions(self) -> pl.DataFrame:
        """Return selected reactions with input and anchor lineage.

        Examples:
            >>> selection.extract_reactions()["ReactionId"].head(1)  # doctest: +SKIP
        """
        return self._extract("reaction", "reaction_id")

    def extract_participants(self) -> pl.DataFrame:
        """Return ordered left/right participants of selected reactions.

        Examples:
            >>> selection.extract_participants().select(  # doctest: +SKIP
            ...     "ReactionId", "Side", "ParticipantId"
            ... )
        """
        return self._extract("reaction_participant", "reaction_id")

    def extract_enzymes(self) -> pl.DataFrame:
        """Return EC links owned by the selected reactions.

        Examples:
            >>> selection.extract_enzymes()["EcNumber"].head(1)  # doctest: +SKIP
        """
        return self._extract("reaction_enzyme", "reaction_id")

    def extract_kos(self) -> pl.DataFrame:
        """Return KO links owned by the selected reactions.

        Examples:
            >>> selection.extract_kos()["KoId"].head(1)  # doctest: +SKIP
        """
        return self._extract("reaction_ko", "reaction_id")

    def extract_modules(self) -> pl.DataFrame:
        """Return module memberships of the selected reactions.

        Examples:
            >>> selection.extract_modules()["ModuleId"].head(1)  # doctest: +SKIP
        """
        return self._extract("reaction_module", "reaction_id")

    def extract_compounds(self) -> pl.DataFrame:
        """Return compound participants and available compound facts.

        Examples:
            >>> selection.extract_compounds()["ParticipantId"].head(1)  # doctest: +SKIP
        """
        participants = self.extract_participants().filter(
            pl.col("ParticipantNamespace") == "kegg_compound"
        )
        self._require("compound")
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            compounds = _query_frame(con, "SELECT * FROM compound")
        return _public_frame(
            participants.join(
                compounds,
                left_on="ParticipantId",
                right_on="compound_id",
                how="left",
            )
        )

    def extract_pathway_memberships(self) -> pl.DataFrame:
        """Return reference-pathway memberships of selected reactions.

        Examples:
            >>> selection.extract_pathway_memberships()[  # doctest: +SKIP
            ...     "PathwayId"
            ... ].head(1)
        """
        return self._extract("reaction_pathway", "reaction_id")

    def extract_cross_references(self) -> pl.DataFrame:
        """Return compound and reaction cross-references in the selection.

        Examples:
            >>> selection.extract_cross_references().select(  # doctest: +SKIP
            ...     "Namespace", "ExternalId"
            ... )
        """
        frames: list[pl.DataFrame] = []
        compounds = self.extract_compounds()
        reactions = self.extract_reactions()
        with duckdb.connect(str(self.publication.path), read_only=True) as con:
            if "compound_cross_reference" in self.publication.tables:
                frames.append(
                    compounds.join(
                        _query_frame(con, "SELECT * FROM compound_cross_reference"),
                        left_on="ParticipantId",
                        right_on="compound_id",
                        how="inner",
                    )
                )
            if "reaction_cross_reference" in self.publication.tables:
                frames.append(
                    reactions.join(
                        _query_frame(con, "SELECT * FROM reaction_cross_reference"),
                        left_on="ReactionId",
                        right_on="reaction_id",
                        how="inner",
                    )
                )
        return (
            _public_frame(pl.concat(frames, how="diagonal_relaxed"))
            if frames
            else pl.DataFrame()
        )

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Return normalized inputs that resolved to no canonical anchor.

        Examples:
            >>> selection.extract_unmatched_ids().select("InputId", "Reason")  # doctest: +SKIP
        """
        matched = {
            (r.get("GroupId"), r["InputId"])
            for r in self.extract_matches().iter_rows(named=True)
        }
        rows: list[dict[str, str | None]] = []
        with duckdb.connect(str(self.publication.path), read_only=True) as connection:
            for group, raw in self.inputs:
                value = _normalize_input(raw, self.namespace)
                if (group, value) in matched:
                    continue
                reason = "not_found"
                if self.namespace == "ec" and "enzyme" in self.publication.tables:
                    status_row = connection.execute(
                        "SELECT status FROM enzyme WHERE ec_number=?", [value]
                    ).fetchone()
                    if status_row is not None and status_row[0] == "deleted":
                        reason = (
                            "not_found"
                            if self.include_obsolete
                            else "obsolete_excluded"
                        )
                    elif status_row is not None and status_row[0] == "transferred":
                        reason = "invalid_canonical_target"
                rows.append({"GroupId": group, "InputId": value, "Reason": reason})
        frame = pl.DataFrame(
            rows,
            schema={"GroupId": pl.String, "InputId": pl.String, "Reason": pl.String},
        )
        return frame if self.is_grouped else frame.drop("GroupId")


def _normalize_input(value: str, namespace: str) -> str:
    value = str(value).strip()
    prefixes = {
        "kegg_compound": "cpd:",
        "kegg_reaction": "rn:",
        "kegg_module": "md:",
        "kegg_pathway": "path:",
    }
    value = value.removeprefix(prefixes.get(namespace, ""))
    if namespace == "chebi":
        value = f"CHEBI:{value.upper().removeprefix('CHEBI:')}"
    if namespace == "rhea":
        value = f"RHEA:{value.upper().removeprefix('RHEA:')}"
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
        raise KEGGMetabolicCapabilityError(
            f"KEGG metabolic namespace {namespace!r} is unavailable; "
            f"available namespaces: {sorted(available)}"
        )


def evaluate_modules(
    publication: MetabolicPublication, ko_ids: Iterable[str]
) -> pl.DataFrame:
    if "module_definition_node" not in publication.tables:
        raise KEGGMetabolicCapabilityError(
            "KEGG metabolic publication lacks module definitions"
        )
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
                "ModuleId": module,
                "RequiredBlockCount": len(blocks),
                "SatisfiedBlockCount": sum(hits),
                "IsComplete": all(hits),
                "MissingBlockIndexes": [i for i, hit in enumerate(hits, 1) if not hit],
            }
        )
    return pl.DataFrame(result).sort("ModuleId")
