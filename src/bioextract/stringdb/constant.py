import re
from typing import Literal

import polars as pl

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
    "input_id": pl.String,
}
SCHEMA_EDGES = {
    "string_id_a": pl.String,
    "string_id_b": pl.String,
    "combined_score": pl.Int64,
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
SCHEMA_GROUP_STRING_MAPPING = {
    "group_id": pl.String,
    "input_id": pl.String,
}
SCHEMA_GROUP_EDGES = {
    "group_id": pl.String,
    "string_id_a": pl.String,
    "string_id_b": pl.String,
    "combined_score": pl.Int64,
}

RE_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")

StringDatabaseVersion = Literal["v12.0"]
StringNamespace = Literal["alias", "string"]
