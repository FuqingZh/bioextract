import re
from dataclasses import fields

from polars._typing import SchemaDict

from .model import (
    AltIdColumnBuffer,
    AncestorColumnBuffer,
    DepthColumnBuffer,
    EdgeColumnBuffer,
    SynonymColumnBuffer,
    TermColumnBuffer,
    XrefColumnBuffer,
)

LARGE_DISTANCE = 1 << 30
SCHEMA_VERSION = "go-obo-tidy-v0.1"
MEDIA_TYPE_OBO = "text/obo"
GO_NAMESPACE_VALUES = (
    "biological_process",
    "molecular_function",
    "cellular_component",
)

RE_GO_ID = re.compile(r"^GO:\d{7}$")
RE_DEFINITION = re.compile(r'^"(?P<text>(?:[^"\\]|\\.)*)"')
RE_SYNONYM = re.compile(
    r'^"(?P<text>(?:[^"\\]|\\.)*)" '
    r"(?P<scope>EXACT|BROAD|NARROW|RELATED)"
    r"(?: (?P<type>[^\[]+?))?"
    r" (?P<dbxref>\[.*\])$"
)

DEDUP_KEYS_BY_FRAME: dict[str, tuple[str, ...]] = {
    "term": ("go_id",),
    "edge": ("child_go_id", "parent_go_id", "relation_type"),
    "synonym": (
        "go_id",
        "synonym_text",
        "synonym_scope",
        "synonym_type_name",
        "dbxref_text",
    ),
    "xref": ("go_id", "xref_text"),
    "alt_id": ("alt_go_id",),
}

ASSET_SPECS: tuple[tuple[str, str, str], ...] = (
    ("term.parquet", "canonical", "term"),
    ("edge.parquet", "canonical", "edge"),
    ("synonym.parquet", "canonical", "synonym"),
    ("xref.parquet", "canonical", "xref"),
    ("alt_id.parquet", "canonical", "alt_id"),
    ("ancestor_all.parquet", "derived", "ancestor_all"),
    ("depth.parquet", "derived", "depth"),
)


SCHEMA_TERM: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(TermColumnBuffer)
}
SCHEMA_EDGE: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(EdgeColumnBuffer)
}
SCHEMA_SYNONYM: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(SynonymColumnBuffer)
}
SCHEMA_XREF: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(XrefColumnBuffer)
}
SCHEMA_ALT_ID: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(AltIdColumnBuffer)
}
SCHEMA_ANCESTOR: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(AncestorColumnBuffer)
}
SCHEMA_DEPTH: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(DepthColumnBuffer)
}
