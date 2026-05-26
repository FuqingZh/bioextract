from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from bioextract._shared import validate_file_size
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .ontology.constant import (
    ASSET_SPECS,
    MEDIA_TYPE_OBO,
    SCHEMA_VERSION,
)
from .ontology.parse import scan_obo_term_records
from .ontology.tidy import build_tidy_frames

__all__ = [
    "GoDb",
    "GoResourceLimits",
    "GoTidyDataset",
]


@dataclass(frozen=True, slots=True)
class GoResourceLimits:
    file_obo_bytes_max: int | None = None


@dataclass(frozen=True, slots=True)
class _GoSnapshot:
    file_obo: Path


GoTidyDataset = TidyDataset


@dataclass(slots=True)
class GoDb:
    """Path-first access to a local Gene Ontology OBO snapshot."""

    snapshot: _GoSnapshot
    limits: GoResourceLimits

    DEFAULT_RESOURCE_LIMITS = GoResourceLimits()

    @classmethod
    def from_obo(
        cls,
        file_obo: os.PathLike[str] | str,
        *,
        limits: GoResourceLimits | None = None,
    ) -> GoDb:
        file_obo = Path(file_obo)
        if not file_obo.exists():
            raise FileNotFoundError(f"GO OBO file not found: {file_obo}")

        limits_resolved = GoResourceLimits() if limits is None else limits
        validate_file_size(
            file_path=file_obo,
            size_max=limits_resolved.file_obo_bytes_max,
            label="GO OBO file",
        )
        return cls(snapshot=_GoSnapshot(file_obo=file_obo), limits=limits_resolved)

    def build_tidy(self) -> GoTidyDataset:
        records = scan_obo_term_records(self.snapshot.file_obo)
        frames = build_tidy_frames(records)
        return GoTidyDataset(
            frames=frames,
            source=TidySource(path=self.snapshot.file_obo, media_type=MEDIA_TYPE_OBO),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"go-ontology-{self.snapshot.file_obo.stem}",
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
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
        )
