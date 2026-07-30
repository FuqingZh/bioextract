from pathlib import Path

import polars as pl

from .constant import (
    DEDUP_KEYS_BY_FRAME,
    SCHEMA_BRITE,
)
from .parse import read_brite
from .write import build_deduped_frame


def build_tidy_frames(file_in: Path) -> dict[str, pl.DataFrame]:
    """Build the canonical tidy frames from one KEGG BRITE hierarchy.

    Args:
        file_in: Path to a KEGG BRITE JSON file.

    Returns:
        A mapping containing the deduplicated ``pathway`` frame. Leafless
        level-three pathways are retained with nullable entry and KO fields.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the JSON hierarchy or a recognized BRITE label is
            malformed.

    Examples:
        Parse the pathway, entry, and KO identifiers from a compact hierarchy:

        >>> frames = build_tidy_frames(Path("data/kegg/tcar00001.json"))
        >>> frames["pathway"].select(
        ...     "pathway_level3_kegg_id", "entry_id", "ko_id"
        ... ).row(0, named=True)
        {'pathway_level3_kegg_id': 'tcar00010', 'entry_id': 'U0034_04525', 'ko_id': 'K00845'}
    """
    brite_buffer = read_brite(file_in)
    df_pathway = build_deduped_frame(
        brite_buffer,
        schema=SCHEMA_BRITE,
        dedup_keys=DEDUP_KEYS_BY_FRAME["pathway"],
    )
    return {"pathway": df_pathway}
