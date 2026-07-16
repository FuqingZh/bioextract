from pathlib import Path

import polars as pl

from .constant import (
    DEDUP_KEYS_BY_FRAME,
    SCHEMA_BRITE,
)
from .parse import read_brite
from .write import build_deduped_frame, write_frame_assets


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
        Build the pathway frame from a compact local hierarchy:

        >>> frames = build_tidy_frames(Path("data/kegg/tcar00001.json"))
        >>> (sorted(frames), frames["pathway"].height)
        (['pathway'], 2)
    """
    brite_buffer = read_brite(file_in)
    df_pathway = build_deduped_frame(
        brite_buffer,
        schema=SCHEMA_BRITE,
        dedup_keys=DEDUP_KEYS_BY_FRAME["pathway"],
    )
    return {"pathway": df_pathway}


def run_tidy_kegg_brite(file_in: Path, dir_out: Path) -> None:
    """Write the legacy manifest-free KEGG BRITE parquet assets.

    Args:
        file_in: Path to a KEGG BRITE JSON file.
        dir_out: Destination directory for ``pathway.parquet``.

    Notes:
        New callers that need source provenance or a manifest should use
        :meth:`bioextract.kegg.KeggDb.write_tidy`. This function remains the
        direct writer for compatibility and intentionally emits no manifest.

    Examples:
        Write the compact BRITE fixture to a relative output directory:

        >>> dir_out = Path("build/kegg-brite")
        >>> run_tidy_kegg_brite(Path("data/kegg/tcar00001.json"), dir_out)
        >>> (dir_out / "pathway.parquet").is_file()
        True
    """
    frames = build_tidy_frames(file_in)
    dir_out.mkdir(parents=True, exist_ok=True)
    write_frame_assets(dir_out=dir_out, frames=frames)
