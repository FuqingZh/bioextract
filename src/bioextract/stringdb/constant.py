import polars as pl
import re
from typing import Literal

DEFAULT_SOURCE_RANK_MAP: dict[str, int] = {
    "UniProt_AC": 1,
    "UniProt_ID": 2,
    "UniProt_GN_Name": 3,
    "UniProt_GN_Synonyms": 4,
    "UniProt_DR_GeneID": 5,
    "KEGG_GENEID": 6,
}

SCHEMA_ALIASES = {
    "v12.0": {
        "#string_protein_id": pl.String,
        "string_protein_id": pl.String,
        "alias": pl.String,
        "source": pl.String,
    },
}
SCHEMA_LINKS = {
    "v12.0": {
        "protein1": pl.String,
        "protein2": pl.String,
        "combined_score": pl.Int64,
    },
}
SCHEMA_PROTEIN_MAP = {
    "InputId": pl.String,
    "StringId": pl.String,
    "MapSource": pl.String,
}
SCHEMA_EDGES = {
    "StringIdA": pl.String,
    "StringIdB": pl.String,
    "Score": pl.Int64,
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
SCHEMA_GROUP_STRING_MAPPING = {
    "GroupId": pl.String,
    "InputId": pl.String,
    "StringId": pl.String,
    "MapSource": pl.String,
}
SCHEMA_GROUP_EDGES = {
    "GroupId": pl.String,
    "StringIdA": pl.String,
    "StringIdB": pl.String,
    "Score": pl.Int64,
}

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")

StringDbVersion = Literal["v12.0"]
