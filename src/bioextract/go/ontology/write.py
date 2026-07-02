import json
from dataclasses import asdict
from pathlib import Path

import polars as pl
from polars._typing import SchemaDict

from .constant import ASSET_SPECS
from .model import FrameColumnBuffer


# #region OutputWriters
def serialize_pretty_json(data: object) -> bytes:
    return (json.dumps(data, indent=2, sort_keys=True) + "\n").encode("utf-8")


def build_deduped_frame(
    data_by_cols: FrameColumnBuffer,
    schema: SchemaDict,
    dedup_keys: tuple[str, ...],
) -> pl.DataFrame:
    return pl.DataFrame(asdict(data_by_cols), schema=schema).unique(
        subset=list(dedup_keys),
        maintain_order=True,
    )


def write_frame_assets(
    frames: dict[str, pl.DataFrame],
    dir_out: Path,
) -> list[dict[str, object]]:
    entries_manifest: list[dict[str, object]] = []
    for path_rel, kind, frame_name in ASSET_SPECS:
        frame = frames[frame_name]
        file_out = dir_out / path_rel
        frame.lazy().sink_parquet(file_out)
        entries_manifest.append(
            {
                "path": path_rel,
                "kind": kind,
                "is_optional": False,
            }
        )

    return entries_manifest


# #endregion
