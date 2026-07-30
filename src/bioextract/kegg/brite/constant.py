from dataclasses import fields

from polars._typing import SchemaDict

from .model import BriteColumnBuffer

SCHEMA_VERSION = "kegg-brite-tidy-v0.1"
MEDIA_TYPE_JSON = "application/json"

DEDUP_KEYS_BY_FRAME: dict[str, tuple[str, ...]] = {
    "pathway": ("pathway_level3_id", "entry_id", "entry_name", "ko_id"),
}

ASSET_SPECS: tuple[tuple[str, str, str], ...] = (
    ("pathway.parquet", "canonical", "pathway"),
)

SCHEMA_BRITE: SchemaDict = {
    field_def.name: field_def.metadata["dtype"]
    for field_def in fields(BriteColumnBuffer)
}
