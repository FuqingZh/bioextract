from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .constant import (
    ASSET_SPECS,
    MEDIA_TYPE_GMT,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    extract_mapping_frame,
    extract_unmapped_input_ids_frame,
    read_gmt_frames,
)

__all__ = [
    "WikiPathwaysDb",
    "WikiPathwaysResourceLimits",
    "WikiPathwaysTidyDataset",
]


@dataclass(frozen=True, slots=True)
class WikiPathwaysResourceLimits:
    file_gmt_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _WikiPathwaysSnapshot:
    file_gmt: Path
    species: str | None = None


WikiPathwaysTidyDataset = TidyDataset


@dataclass(slots=True)
class WikiPathwaysDb:
    """Path-first access to a local WikiPathways GMT snapshot.

    `WikiPathwaysDb` is the public entrypoint for extracting pathway gene sets
    from local WikiPathways GMT files. GMT rows are interpreted as pathway
    metadata followed by NCBI Entrez Gene IDs; the class does not perform
    identifier conversion or enrichment statistics.

    Construct instances with :meth:`from_gmt`, then either extract whole
    resource frames or create single/grouped Entrez ID selections through
    :meth:`select_ids` and :meth:`select_groups`.
    """

    snapshot: _WikiPathwaysSnapshot
    limits: WikiPathwaysResourceLimits = field(
        default_factory=WikiPathwaysResourceLimits
    )
    _frames: dict[str, pl.LazyFrame] | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_RESOURCE_LIMITS = WikiPathwaysResourceLimits()

    @classmethod
    def from_gmt(
        cls,
        file_gmt: os.PathLike[str] | str,
        *,
        species: str | None = None,
        limits: WikiPathwaysResourceLimits | None = None,
    ) -> WikiPathwaysDb:
        """Create a dataset handle from a local WikiPathways GMT file.

        Args:
            file_gmt: Path to a local WikiPathways GMT file.
            species: Optional species display name used as an exact metadata
                filter after parsing.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used.

        Returns:
            A dataset handle that can build pathway, term2gene, and term2name
            frames.

        Raises:
            FileNotFoundError: If the GMT file does not exist.
            ValueError: If the species string is empty after normalization or
                a configured file-size limit is exceeded.
        """
        file_gmt = Path(file_gmt)
        if not file_gmt.exists():
            raise FileNotFoundError(f"WikiPathways GMT file not found: {file_gmt}")
        limits_resolved = WikiPathwaysResourceLimits() if limits is None else limits
        validate_file_size(
            file_path=file_gmt,
            size_max=limits_resolved.file_gmt_bytes_max,
            label="WikiPathways GMT file",
        )
        species_normalized = None if species is None else str(species).strip()
        if species is not None and not species_normalized:
            raise ValueError(
                "WikiPathways species must be non-empty after normalization"
            )
        return cls(
            snapshot=_WikiPathwaysSnapshot(
                file_gmt=file_gmt,
                species=species_normalized,
            ),
            limits=limits_resolved,
        )

    def select_ids(self, ids: Iterable[str]) -> WikiPathwaysSelection:
        """Create a single-query selection from Entrez Gene IDs.

        Args:
            ids: Input Entrez Gene IDs. Values are normalized with the shared
                input-ID normalizer.

        Returns:
            A selection that can extract pathway mappings and unmapped IDs.

        Raises:
            ValueError: If the normalized input-ID count exceeds the configured
                limit.
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return WikiPathwaysSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
    ) -> WikiPathwaysSelection:
        """Create a grouped selection from multiple Entrez Gene ID sets.

        Args:
            group_to_ids: Mapping from group label to input Entrez Gene IDs.

        Returns:
            A grouped selection that carries `GroupId` through outputs.

        Raises:
            ValueError: If group or input-ID limits are exceeded, or if group
                IDs are invalid after normalization.
        """
        grp_in_frames = create_group_input_frames(
            group_to_ids,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        validate_count_limit(
            count=grp_in_frames.df_groups.height,
            limit_max=self.limits.num_groups_max,
            label="Group count",
        )
        validate_count_limit(
            count=grp_in_frames.df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return WikiPathwaysSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
        )

    def extract_pathway(self) -> pl.DataFrame:
        """Extract WikiPathways pathway metadata with gene counts."""
        return self._lazy_frame("pathway").collect()

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract a WikiPathways pathway-to-Entrez table."""
        return self._lazy_frame("term2gene").collect()

    def extract_term2name(self) -> pl.DataFrame:
        """Extract WikiPathways pathway display metadata for enrichment callers."""
        return self._lazy_frame("term2name").collect()

    def build_tidy(self) -> WikiPathwaysTidyDataset:
        """Build the in-memory WikiPathways tidy dataset.

        Returns:
            A `TidyDataset` with `pathway`, `term2gene`, and `term2name` frames.
        """
        return WikiPathwaysTidyDataset(
            frames={
                "pathway": self._lazy_frame("pathway"),
                "term2gene": self._lazy_frame("term2gene"),
                "term2name": self._lazy_frame("term2name"),
            },
            source=TidySource(path=self.snapshot.file_gmt, media_type=MEDIA_TYPE_GMT),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"wikipathways-gmt-{self.snapshot.file_gmt.stem}",
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
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the WikiPathways tidy dataset as flat parquet files.

        Args:
            dir_out: Output directory for parquet assets.
            should_write_manifest: Whether to write `manifest.json`.
            should_hash_assets: Whether to calculate asset checksums in the
                manifest.

        Returns:
            A write report with asset paths and optional manifest content.
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def _lazy_frame(self, frame_name: str) -> pl.LazyFrame:
        if self._frames is None:
            frames = read_gmt_frames(self.snapshot.file_gmt)
            if self.snapshot.species is not None:
                lf_pathway = _filter_species_frame(
                    frames["pathway"],
                    self.snapshot.species,
                )
                lf_pathway_ids = lf_pathway.select("WikiPathwaysId").unique()
                frames = {
                    "pathway": lf_pathway,
                    "term2gene": frames["term2gene"]
                    .join(
                        lf_pathway_ids,
                        on="WikiPathwaysId",
                        how="inner",
                    )
                    .sort("WikiPathwaysId", "GeneId"),
                    "term2name": _filter_species_frame(
                        frames["term2name"],
                        self.snapshot.species,
                    ),
                }
            self._frames = frames
        return self._frames[frame_name]


@dataclass(slots=True)
class WikiPathwaysSelection:
    """Selection handle for single and grouped WikiPathways queries.

    Selections are created by :meth:`WikiPathwaysDb.select_ids` or
    :meth:`WikiPathwaysDb.select_groups`. Single selections return tables keyed
    by `InputId`; grouped selections prepend `GroupId`.
    """

    dataset: WikiPathwaysDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _lf_mapping: pl.LazyFrame | None = field(default=None, repr=False)
    _lf_unmapped: pl.LazyFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs."""
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the group ID column when this selection is grouped."""
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract selected Entrez-to-WikiPathways pathway mappings."""
        return self._lazy_mapping().collect()

    def _lazy_mapping(self) -> pl.LazyFrame:
        if self._lf_mapping is None:
            self._lf_mapping = extract_mapping_frame(
                self.dataset._lazy_frame("pathway"),
                self.dataset._lazy_frame("term2gene"),
                self._df_input_ids,
                cols_group_id=self._col_group_id,
            )
        return self._lf_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs that did not map to WikiPathways."""
        if self._lf_unmapped is None:
            self._lf_unmapped = extract_unmapped_input_ids_frame(
                self._df_input_ids,
                self._lazy_mapping(),
                cols_group_id=self._col_group_id,
            )
        return self._lf_unmapped.collect()


def _filter_species_frame(lf: pl.LazyFrame, species: str) -> pl.LazyFrame:
    if "Species" not in lf.collect_schema().names():
        return lf
    return lf.filter(pl.col("Species") == species)
