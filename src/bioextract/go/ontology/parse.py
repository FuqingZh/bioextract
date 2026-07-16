from collections import defaultdict
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

from .constant import RE_DEFINITION, RE_GO_ID, RE_SYNONYM
from .model import (
    ParentEdge,
    SubsetDefinitionRecord,
    SynonymRecord,
    TermRecord,
    XrefParts,
)


# #region OboParsing
def extract_inline_value(value: str) -> str:
    return value.split(" !", 1)[0].strip()


def normalize_whitespace(text: str) -> str:
    return " ".join(text.strip().split())


def validate_go_id(go_id: str) -> None:
    if not RE_GO_ID.match(go_id):
        raise ValueError(f"Invalid GO identifier: {go_id!r}")


def parse_bool_text(value: str | None, *, default: bool = False) -> bool:
    if value is None:
        return default
    return value.strip().lower() == "true"


def parse_definition_text(raw_definition: str | None) -> str | None:
    if raw_definition is None:
        return None
    match = RE_DEFINITION.match(raw_definition.strip())
    if match is None:
        return None
    return match.group("text")


def parse_synonym(synonym: str) -> SynonymRecord | None:
    synonym_trimmed = normalize_whitespace(synonym)
    match = RE_SYNONYM.match(synonym_trimmed)
    if match is None:
        return None

    synonym_type_name = (
        None if (value := match.group("type")) is None else value.strip() or None
    )
    dbxref_text = (
        block[1:-1].strip()
        if (block := match.group("dbxref")) not in {"[]", None}
        else None
    )
    return SynonymRecord(
        synonym_text=match.group("text"),
        synonym_scope=match.group("scope"),
        synonym_type_name=synonym_type_name,
        dbxref_text=dbxref_text,
    )


def parse_xref_lossless(xref: str) -> XrefParts:
    """Split only unambiguous ``DB:ID`` xrefs into structured fields.

    Whitespace, quotes, missing components, or a missing colon make the xref
    unsafe to normalize. Those values return two ``None`` fields so the tidy
    layer can keep the original xref text without inventing structure.
    """
    is_empty = not xref
    is_whitespace_present = any(char.isspace() for char in xref)
    is_quote_present = '"' in xref
    is_colon_present = ":" in xref
    if is_empty or is_whitespace_present or is_quote_present or not is_colon_present:
        return XrefParts(None, None)

    xref_db, xref_id = xref.split(":", 1)
    if not xref_db or not xref_id:
        return XrefParts(None, None)

    return XrefParts(xref_db, xref_id)


def parse_subset_definition(raw_subsetdef: str) -> SubsetDefinitionRecord | None:
    subset_id, _, raw_subset_name = raw_subsetdef.partition(" ")
    subset_id = subset_id.strip()
    subset_name = raw_subset_name.strip()
    if not subset_id or not subset_name:
        return None
    if subset_name.startswith('"') and subset_name.endswith('"'):
        subset_name = subset_name[1:-1]
    return SubsetDefinitionRecord(subset_id=subset_id, subset_name=subset_name)


def parse_parent_edges(
    block: dict[str, list[str]],
    child_go_id: str,
) -> list[ParentEdge]:
    """Parse ``is_a`` and generic OBO relationship clauses as parent edges.

    Malformed generic relationship payloads are skipped, while syntactically
    present parent IDs are validated as GO identifiers. ``source_clause``
    remains distinct from ``relation_type`` so downstream graph policy can
    distinguish OBO syntax from biological relation semantics.
    """
    parents: list[ParentEdge] = []
    for raw_parent in block.get("is_a", []):
        parent_go_id = extract_inline_value(raw_parent)
        validate_go_id(parent_go_id)
        parents.append(
            ParentEdge(
                child_go_id=child_go_id,
                parent_go_id=parent_go_id,
                relation_type="is_a",
                source_clause="is_a",
            )
        )

    # OBO relationship lines are "<relation_type> <parent_go_id>" before the inline comment.
    for raw_relationship in block.get("relationship", []):
        relationship_payload = extract_inline_value(raw_relationship)
        relationship_components = relationship_payload.split(maxsplit=1)
        if len(relationship_components) != 2:
            continue

        relation_type, parent_go_id = relationship_components
        validate_go_id(parent_go_id)
        parents.append(
            ParentEdge(
                child_go_id=child_go_id,
                parent_go_id=parent_go_id,
                relation_type=relation_type,
                source_clause="relationship",
            )
        )

    return parents


def parse_term_record(block: dict[str, list[str]]) -> TermRecord | None:
    go_id = extract_inline_value(values[-1]) if (values := block.get("id")) else None
    term_name = (
        extract_inline_value(values[-1]) if (values := block.get("name")) else None
    )
    namespace = (
        extract_inline_value(values[-1]) if (values := block.get("namespace")) else None
    )
    if go_id is None or term_name is None or namespace is None:
        return None

    validate_go_id(go_id)
    definition = parse_definition_text(
        extract_inline_value(values[-1]) if (values := block.get("def")) else None
    )
    comment = (
        extract_inline_value(values[-1]) if (values := block.get("comment")) else None
    )
    is_obsolete = parse_bool_text(
        extract_inline_value(values[-1])
        if (values := block.get("is_obsolete"))
        else None
    )

    return TermRecord(
        go_id=go_id,
        term_name=term_name,
        namespace=namespace,
        definition=definition,
        is_obsolete=is_obsolete,
        comment=comment,
        alt_ids=[extract_inline_value(value) for value in block.get("alt_id", [])],
        subsets=[extract_inline_value(value) for value in block.get("subset", [])],
        xrefs=[extract_inline_value(value) for value in block.get("xref", [])],
        synonyms=[value.strip() for value in block.get("synonym", [])],
        parents=parse_parent_edges(block=block, child_go_id=go_id),
    )


def scan_obo_term_records(file_in: Path) -> Iterator[TermRecord]:
    """Stream complete ``[Term]`` records from a GO OBO snapshot.

    Non-term stanzas and term blocks missing ID, name, or namespace are
    ignored. EOF is treated as a stanza boundary, so a final term does not
    require a trailing blank line.
    """
    current_stanza_kind: str | None = None
    block: dict[str, list[str]] = defaultdict(list)

    with file_in.open("r", encoding="utf-8") as handle:
        for raw_line in chain(handle, ["\n"]):  # Treat EOF as a final blank line.
            line_stripped = raw_line.rstrip("\n").strip()
            is_blank_line = not line_stripped
            is_stanza_header = line_stripped.startswith("[") and line_stripped.endswith(
                "]"
            )
            should_commit_block = is_blank_line or is_stanza_header
            is_term_stanza = current_stanza_kind == "Term"
            is_field_line = ": " in line_stripped

            if should_commit_block and is_term_stanza and block:
                record = parse_term_record(block)
                if record is not None:
                    yield record

            if is_blank_line:
                current_stanza_kind = None
                block = defaultdict(list)
                continue

            if is_stanza_header:
                current_stanza_kind = line_stripped[1:-1]
                block = defaultdict(list)
                continue

            if not is_term_stanza or not is_field_line:
                continue

            key, value = line_stripped.split(": ", 1)
            block[key].append(value.strip())


def read_obo_term_records(file_in: Path) -> list[TermRecord]:
    return list(scan_obo_term_records(file_in))


def read_obo_subset_definitions(file_in: Path) -> list[SubsetDefinitionRecord]:
    definitions: list[SubsetDefinitionRecord] = []
    with file_in.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line_stripped = raw_line.rstrip("\n").strip()
            if line_stripped.startswith("["):
                break
            if not line_stripped.startswith("subsetdef: "):
                continue
            definition = parse_subset_definition(
                line_stripped.removeprefix("subsetdef: ")
            )
            if definition is not None:
                definitions.append(definition)
    return definitions


# #endregion
