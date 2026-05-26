from pathlib import Path

import polars as pl

from .constant import (
    DEDUP_KEYS_BY_FRAME,
    SCHEMA_BRITE,
)
from .parse import read_brite
from .write import build_deduped_frame, write_frame_assets


def build_tidy_frames(file_in: Path) -> dict[str, pl.DataFrame]:
    brite_buffer = read_brite(file_in)
    df_pathway = build_deduped_frame(
        brite_buffer,
        schema=SCHEMA_BRITE,
        dedup_keys=DEDUP_KEYS_BY_FRAME["pathway"],
    )
    return {"pathway": df_pathway}


def run_tidy_kegg_brite(file_in: Path, dir_out: Path) -> None:
    frames = build_tidy_frames(file_in)
    dir_out.mkdir(parents=True, exist_ok=True)
    write_frame_assets(dir_out=dir_out, frames=frames)
