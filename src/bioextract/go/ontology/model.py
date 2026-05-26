from __future__ import annotations

from dataclasses import dataclass, field

import polars as pl


@dataclass(slots=True)
class XrefParts:
    xref_db: str | None
    xref_id: str | None


@dataclass(slots=True)
class ParentEdge:
    child_go_id: str
    parent_go_id: str
    relation_type: str
    source_clause: str | None = None


@dataclass(slots=True)
class SynonymRecord:
    synonym_text: str
    synonym_scope: str
    synonym_type_name: str | None
    dbxref_text: str | None


@dataclass(slots=True)
class TermRecord:
    go_id: str
    term_name: str
    namespace: str
    definition: str | None
    is_obsolete: bool
    comment: str | None
    alt_ids: list[str] = field(default_factory=lambda: [])
    xrefs: list[str] = field(default_factory=lambda: [])
    synonyms: list[str] = field(default_factory=lambda: [])
    parents: list[ParentEdge] = field(default_factory=lambda: [])


@dataclass(slots=True)
class TermColumnBuffer:
    go_id: list[str] = field(default_factory=lambda: [], metadata={"dtype": pl.String})
    term_name: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    namespace: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    definition: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    is_obsolete: list[bool] = field(
        default_factory=lambda: [], metadata={"dtype": pl.Boolean}
    )
    comment: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def append_record(self, record: TermRecord) -> None:
        self.go_id.append(record.go_id)
        self.term_name.append(record.term_name)
        self.namespace.append(record.namespace)
        self.definition.append(record.definition)
        self.is_obsolete.append(record.is_obsolete)
        self.comment.append(record.comment)


@dataclass(slots=True)
class EdgeColumnBuffer:
    child_go_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    parent_go_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    relation_type: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    source_clause: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def append_edge(self, edge: ParentEdge) -> None:
        self.child_go_id.append(edge.child_go_id)
        self.parent_go_id.append(edge.parent_go_id)
        self.relation_type.append(edge.relation_type)
        self.source_clause.append(edge.source_clause)


@dataclass(slots=True)
class AltIdColumnBuffer:
    alt_go_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    primary_go_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def extend(self, other: AltIdColumnBuffer) -> None:
        self.alt_go_id.extend(other.alt_go_id)
        self.primary_go_id.extend(other.primary_go_id)


@dataclass(slots=True)
class SynonymColumnBuffer:
    go_id: list[str] = field(default_factory=lambda: [], metadata={"dtype": pl.String})
    synonym_text: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    synonym_scope: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    synonym_type_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    dbxref_text: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def extend(self, other: SynonymColumnBuffer) -> None:
        self.go_id.extend(other.go_id)
        self.synonym_text.extend(other.synonym_text)
        self.synonym_scope.extend(other.synonym_scope)
        self.synonym_type_name.extend(other.synonym_type_name)
        self.dbxref_text.extend(other.dbxref_text)


@dataclass(slots=True)
class XrefColumnBuffer:
    go_id: list[str] = field(default_factory=lambda: [], metadata={"dtype": pl.String})
    xref_text: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    xref_db: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    xref_id: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def extend(self, other: XrefColumnBuffer) -> None:
        self.go_id.extend(other.go_id)
        self.xref_text.extend(other.xref_text)
        self.xref_db.extend(other.xref_db)
        self.xref_id.extend(other.xref_id)


@dataclass(slots=True)
class AncestorColumnBuffer:
    go_id: list[str] = field(default_factory=lambda: [], metadata={"dtype": pl.String})
    ancestor_go_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    min_distance: list[int] = field(
        default_factory=lambda: [], metadata={"dtype": pl.Int32}
    )


@dataclass(slots=True)
class DepthColumnBuffer:
    go_id: list[str] = field(default_factory=lambda: [], metadata={"dtype": pl.String})
    namespace: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    min_depth_from_root: list[int] = field(
        default_factory=lambda: [], metadata={"dtype": pl.Int32}
    )
    max_depth_from_root: list[int] = field(
        default_factory=lambda: [], metadata={"dtype": pl.Int32}
    )


FrameColumnBuffer = (
    EdgeColumnBuffer
    | TermColumnBuffer
    | SynonymColumnBuffer
    | XrefColumnBuffer
    | AltIdColumnBuffer
    | AncestorColumnBuffer
    | DepthColumnBuffer
)
