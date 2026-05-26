from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from bioextract._shared import validate_file_size
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .brite.constant import (
    ASSET_SPECS,
    MEDIA_TYPE_JSON,
    SCHEMA_VERSION,
)
from .brite.tidy import build_tidy_frames

__all__ = [
    "KeggDb",
    "KeggResourceLimits",
    "KeggTidyDataset",
]


@dataclass(frozen=True, slots=True)
class KeggResourceLimits:
    file_brite_json_bytes_max: int | None = None


@dataclass(frozen=True, slots=True)
class _KeggSnapshot:
    file_brite_json: Path


KeggTidyDataset = TidyDataset


@dataclass(slots=True)
class KeggDb:
    """Path-first access to local KEGG resource snapshots.

    `KeggDb` is the public entrypoint for extracting tidy KEGG BRITE pathway
    tables from a local JSON snapshot. It stores only the source path and
    resource limits until `build_tidy()` or `write_tidy()` is called.
    """

    snapshot: _KeggSnapshot
    limits: KeggResourceLimits

    DEFAULT_RESOURCE_LIMITS = KeggResourceLimits()

    @classmethod
    def from_brite_json(
        cls,
        file_brite_json: os.PathLike[str] | str,
        *,
        limits: KeggResourceLimits | None = None,
    ) -> KeggDb:
        """Create a dataset handle from a local KEGG BRITE JSON file.

        Args:
            file_brite_json: Path to a local KEGG BRITE JSON resource.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used.

        Returns:
            A dataset handle that can build tidy KEGG pathway frames.

        Raises:
            FileNotFoundError: If the BRITE JSON file does not exist.
            ValueError: If the configured file-size limit is exceeded.
        """
        file_brite_json = Path(file_brite_json)
        if not file_brite_json.exists():
            raise FileNotFoundError(f"KEGG BRITE JSON file not found: {file_brite_json}")

        limits_resolved = KeggResourceLimits() if limits is None else limits
        validate_file_size(
            file_path=file_brite_json,
            size_max=limits_resolved.file_brite_json_bytes_max,
            label="KEGG BRITE JSON file",
        )
        return cls(
            snapshot=_KeggSnapshot(file_brite_json=file_brite_json),
            limits=limits_resolved,
        )

    def build_tidy(self) -> KeggTidyDataset:
        """Build the in-memory KEGG BRITE tidy dataset.

        Returns:
            A `TidyDataset` containing the flat `pathway` frame.
        """
        frames = build_tidy_frames(self.snapshot.file_brite_json)
        return KeggTidyDataset(
            frames=frames,
            source=TidySource(
                path=self.snapshot.file_brite_json,
                media_type=MEDIA_TYPE_JSON,
            ),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"kegg-brite-{self.snapshot.file_brite_json.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
    ) -> TidyWriteReport:
        """Write the KEGG tidy dataset as flat parquet files.

        Args:
            dir_out: Output directory for parquet assets.
            should_write_manifest: Whether to write `manifest.json`.

        Returns:
            A write report with asset paths and optional manifest content.
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
        )
