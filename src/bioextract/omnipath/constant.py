import re
from typing import Literal

import polars as pl
from polars._typing import SchemaDict

SCHEMA_ENZSUB_RAW = {
    "enzyme": pl.String,
    "substrate": pl.String,
    "residue_type": pl.String,
    "residue_offset": pl.String,
    "modification": pl.String,
}

SCHEMA_INTERACTIONS_RAW = {
    "source": pl.String,
    "target": pl.String,
    "is_directed": pl.String,
    "is_stimulation": pl.String,
    "is_inhibition": pl.String,
}

SCHEMA_UNMAPPED = {
    "input_id": pl.String,
}

SCHEMA_GROUPS = {
    "group_id": pl.String,
}

SCHEMA_GROUP_INPUT_IDS = {
    "group_id": pl.String,
    "input_id": pl.String,
}

SCHEMA_ENZSUB = {
    "enzyme": pl.String,
    "substrate": pl.String,
    "residue_type": pl.String,
    "residue_offset": pl.String,
    "modification": pl.String,
    "target_site": pl.String,
}

SCHEMA_GROUP_ENZSUB = {
    "group_id": pl.String,
    **SCHEMA_ENZSUB,
}

SCHEMA_INTERACTIONS: SchemaDict = {
    "source": pl.String,
    "target": pl.String,
    "is_directed": pl.Boolean,
    "is_stimulation": pl.Boolean,
    "is_inhibition": pl.Boolean,
}

SCHEMA_GROUP_INTERACTIONS = {
    "group_id": pl.String,
    **SCHEMA_INTERACTIONS,
}

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")

OmniPathResourceName = Literal["enzsub", "interactions"]
OmniPathNamespace = Literal["protein"]
