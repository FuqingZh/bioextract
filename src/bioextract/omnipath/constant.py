import re
from typing import Literal

import polars as pl

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
    "InputId": pl.String,
}

SCHEMA_GROUPS = {
    "GroupId": pl.String,
}

SCHEMA_GROUP_INPUT_IDS = {
    "GroupId": pl.String,
    "InputId": pl.String,
}

SCHEMA_ENZSUB = {
    "SourceId": pl.String,
    "TargetId": pl.String,
    "TargetSite": pl.String,
    "Modification": pl.String,
}

SCHEMA_GROUP_ENZSUB = {
    "GroupId": pl.String,
    "SourceId": pl.String,
    "TargetId": pl.String,
    "TargetSite": pl.String,
    "Modification": pl.String,
}

SCHEMA_INTERACTIONS = {
    "SourceId": pl.String,
    "TargetId": pl.String,
    "IsDirected": pl.Boolean,
    "IsStimulation": pl.Boolean,
    "IsInhibition": pl.Boolean,
}

SCHEMA_GROUP_INTERACTIONS = {
    "GroupId": pl.String,
    "SourceId": pl.String,
    "TargetId": pl.String,
    "IsDirected": pl.Boolean,
    "IsStimulation": pl.Boolean,
    "IsInhibition": pl.Boolean,
}

COLS_RENAMED_ENZSUB = {
    "enzyme": "SourceId",
    "substrate": "TargetId",
    "residue_type": "_ResidueType",
    "residue_offset": "_ResidueOffset",
    "modification": "Modification",
}

COLS_RENAMED_INTERACTIONS = {
    "source": "SourceId",
    "target": "TargetId",
    "is_directed": "IsDirected",
    "is_stimulation": "IsStimulation",
    "is_inhibition": "IsInhibition",
}

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")

OmniPathResourceName = Literal["enzsub", "interactions"]
