from pathlib import Path
import json

from .model import (
    BriteColumnBuffer,
    BriteRecord,
    PathwayLeafRecord,
    PathwayLevelRecord,
)


def parse_category(raw: str) -> PathwayLevelRecord:
    parts = raw.strip().split(maxsplit=1)
    if len(parts) == 2 and parts[0].isdigit():
        return PathwayLevelRecord(id=parts[0], name=parts[1])
    raise ValueError(f"Invalid KEGG BRITE category node: {raw!r}")


def parse_pathway_level3(raw: str) -> PathwayLevelRecord:
    text = raw.strip()
    if text.endswith("]") and " [" in text:
        pathway_description, level3_payload = text[:-1].rsplit(" [", 1)
        pathway_parts = pathway_description.split(maxsplit=1)
        if (
            len(pathway_parts) != 2
            or not pathway_parts[0].isdigit()
            or ":" not in level3_payload
        ):
            raise ValueError(f"Invalid KEGG BRITE level-3 pathway node: {raw!r}")

        _, pathway_name = pathway_parts
        _, pathway_id = level3_payload.split(":", 1)
        if not pathway_id:
            raise ValueError(f"Invalid KEGG BRITE level-3 pathway node: {raw!r}")
        return PathwayLevelRecord(
            id=pathway_parts[0],
            name=pathway_name,
            kegg_id=pathway_id,
        )

    pathway_parts = text.split(maxsplit=1)
    if len(pathway_parts) != 2 or not pathway_parts[0].isdigit():
        raise ValueError(f"Invalid KEGG BRITE level-3 pathway node: {raw!r}")

    pathway_id, pathway_name = pathway_parts
    return PathwayLevelRecord(id=pathway_id, name=pathway_name)


def parse_entry_and_ko(raw: str) -> PathwayLeafRecord:
    parts = raw.split("\t", 1)
    if len(parts) == 1:
        entry_tokens = parts[0].strip().split(maxsplit=1)
        if not entry_tokens:
            raise ValueError(f"Invalid KEGG BRITE entry node: {raw!r}")
        entry_id = entry_tokens[0]
        entry_name = entry_tokens[1] if len(entry_tokens) == 2 else None
        return PathwayLeafRecord(
            entry=PathwayLevelRecord(id=entry_id, name=entry_name),
            ko=None,
        )

    if len(parts) != 2:
        raise ValueError(f"Invalid KEGG BRITE entry node: {raw!r}")

    entry_part, ko_part = (part.strip() for part in parts)
    entry_tokens = entry_part.split(maxsplit=1)
    ko_tokens = ko_part.split(maxsplit=1)
    if not entry_tokens or not ko_tokens:
        raise ValueError(f"Invalid KEGG BRITE entry node: {raw!r}")

    entry_id = entry_tokens[0]
    entry_name = entry_tokens[1] if len(entry_tokens) == 2 else None
    ko_id = ko_tokens[0]
    ko_name = ko_tokens[1] if len(ko_tokens) == 2 else None

    return PathwayLeafRecord(
        entry=PathwayLevelRecord(id=entry_id, name=entry_name),
        ko=PathwayLevelRecord(id=ko_id, name=ko_name),
    )


def read_brite(file_in: Path) -> BriteColumnBuffer:
    data = json.loads(file_in.read_text(encoding="utf-8"))
    records = BriteColumnBuffer()

    for level1_node in data.get("children", []):
        pathway_level1 = parse_category(level1_node["name"])
        for level2_node in level1_node.get("children", []):
            pathway_level2 = parse_category(level2_node["name"])
            for level3_node in level2_node.get("children", []):
                pathway_level3 = parse_pathway_level3(level3_node["name"])
                leaf_nodes = level3_node.get("children", [])
                if not leaf_nodes:
                    records.append_record(
                        BriteRecord(
                            pathway_level1_id=pathway_level1.id,
                            pathway_level1_name=pathway_level1.name,
                            pathway_level2_id=pathway_level2.id,
                            pathway_level2_name=pathway_level2.name,
                            pathway_level3_id=pathway_level3.id,
                            pathway_level3_kegg_id=pathway_level3.kegg_id,
                            pathway_level3_name=pathway_level3.name,
                            entry_id=None,
                            entry_name=None,
                            ko_id=None,
                            ko_name=None,
                        )
                    )
                    continue

                for leaf_node in leaf_nodes:
                    leaf_record = parse_entry_and_ko(leaf_node["name"])
                    records.append_record(
                        BriteRecord(
                            pathway_level1_id=pathway_level1.id,
                            pathway_level1_name=pathway_level1.name,
                            pathway_level2_id=pathway_level2.id,
                            pathway_level2_name=pathway_level2.name,
                            pathway_level3_id=pathway_level3.id,
                            pathway_level3_kegg_id=pathway_level3.kegg_id,
                            pathway_level3_name=pathway_level3.name,
                            entry_id=leaf_record.entry.id,
                            entry_name=leaf_record.entry.name,
                            ko_id=leaf_record.ko.id if leaf_record.ko else None,
                            ko_name=leaf_record.ko.name if leaf_record.ko else None,
                        )
                    )

    return records
