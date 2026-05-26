from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from bioextract._shared import validate_file_size
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .ontology.constant import (
    ASSET_SPECS,
    MEDIA_TYPE_OBO,
    SCHEMA_VERSION,
)
from .ontology.parse import scan_obo_term_records
from .ontology.tidy import build_tidy_frames, extract_subcell_frame

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
    """Path-first access to a local Gene Ontology OBO snapshot.

    `GoDb` is the public entrypoint for extracting tidy ontology tables from a
    local GO OBO file. It keeps the raw file path and resource limits, then
    builds materialized Polars frames only when tidy or convenience exports are
    requested.

    The default tidy output is a flat ontology snapshot with canonical term and
    edge tables plus derived graph tables. `extract_subcell()` is a convenience
    view over cellular component terms for subcellular-location workflows.
    """

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
        """Create a dataset handle from a local GO OBO file.

        Args:
            file_obo: Path to a local Gene Ontology OBO file.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used.

        Returns:
            A dataset handle that can build tidy ontology frames and subcellular
            component exports.

        Raises:
            FileNotFoundError: If the OBO file does not exist.
            ValueError: If the configured file-size limit is exceeded.
        """
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
        """Build the in-memory GO tidy dataset.

        Returns:
            A `TidyDataset` with `term`, `edge`, `synonym`, `xref`, `alt_id`,
            `ancestor_all`, and `depth` frames.
        """
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
        """Write the GO tidy dataset as flat parquet files.

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

    def extract_subcell(self, *, include_obsolete: bool = False) -> pl.DataFrame:
        """Extract non-obsolete cellular component terms as a subcell table.

        Args:
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            A DataFrame with GO ID, subcell name, definition, and depth columns.
        """
        return extract_subcell_frame(
            self.build_tidy().frames,
            include_obsolete=include_obsolete,
        )

    def write_subcell(
        self,
        file_out: os.PathLike[str] | str,
        *,
        include_obsolete: bool = False,
    ) -> Path:
        """Write the cellular component subcell table as a parquet file.

        Args:
            file_out: Output parquet path.
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            The output path that was written.
        """
        file_out = Path(file_out)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        self.extract_subcell(include_obsolete=include_obsolete).write_parquet(file_out)
        return file_out
