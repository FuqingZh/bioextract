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
    """Optional fail-fast limits for WikiPathways resources and selections.

    The GMT file limit is measured in bytes and checked when a snapshot handle
    is created. Count limits apply after input IDs and group names are
    normalized; ``None`` disables the corresponding limit.

    Examples:
        Limit one normalized selection to 500 Entrez IDs:

        >>> limits = WikiPathwaysResourceLimits(num_input_ids_max=500)
        >>> limits.num_input_ids_max
        500
    """

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

    Examples:
        Extract pathway mappings while preserving unmapped Entrez IDs:

        >>> db = WikiPathwaysDb.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
        ...     species="Homo sapiens",
        ... )
        >>> db.select_ids(["2687", "MISSING"]).extract_mapping().height
        1
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
            limits: Optional resource policy. When omitted, no finite size or
                selection limits are imposed.

        Returns:
            A dataset handle that can build pathway, term2gene, and term2name
            frames.

        Raises:
            FileNotFoundError: If the GMT file does not exist.
            ValueError: If the species string is empty after normalization or
                a configured file-size limit is exceeded.

        Examples:
            Open a compact human WikiPathways fixture:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_pathway().height
            2
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

        Examples:
            Create an ungrouped Entrez ID selection:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> db.select_ids(["2687"]).is_grouped
            False
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

        Examples:
            Create a grouped Entrez ID selection:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> db.select_groups({"case": ["2687"]}).is_grouped
            True
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
        """Extract one metadata row per pathway in the current species scope.

        ``GeneCount`` counts distinct, non-empty Entrez Gene IDs in the GMT row.

        Examples:
            List pathway IDs in the compact human fixture:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_pathway()["WikiPathwaysId"].to_list()
            ['WP100', 'WP106']
        """
        return self.lazy_frame("pathway").collect()

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract distinct WikiPathways-pathway-to-Entrez enrichment pairs.

        Examples:
            Inspect the compact human enrichment mapping shape:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_term2gene().shape
            (4, 2)
        """
        return self.lazy_frame("term2gene").collect()

    def extract_term2name(self) -> pl.DataFrame:
        """Extract one pathway display-metadata row per WikiPathways ID.

        Examples:
            Inspect the enrichment label schema:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> db.extract_term2name().columns
            ['WikiPathwaysId', 'PathwayName', 'Species', 'Collection', 'Version', 'Url']
        """
        return self.lazy_frame("term2name").collect()

    def build_tidy(self) -> WikiPathwaysTidyDataset:
        """Build the lazy WikiPathways tidy dataset for the current scope.

        Returns:
            A `TidyDataset` with species-consistent ``pathway``, ``term2gene``,
            and ``term2name`` frames.

        Examples:
            Build the three declared WikiPathways frames:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> sorted(db.build_tidy().frames)
            ['pathway', 'term2gene', 'term2name']
        """
        return WikiPathwaysTidyDataset(
            frames={
                "pathway": self.lazy_frame("pathway"),
                "term2gene": self.lazy_frame("term2gene"),
                "term2name": self.lazy_frame("term2name"),
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
                manifest. This has no effect unless a manifest is requested.

        Returns:
            A write report with asset paths and optional manifest content.

        Examples:
            Write the three declared WikiPathways assets:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> report = db.write_tidy("build/wikipathways")
            >>> [asset.path for asset in report.assets]
            ['pathway.parquet', 'term2gene.parquet', 'term2name.parquet']
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def lazy_frame(self, frame_name: str) -> pl.LazyFrame:
        """Return a cached lazy tidy frame in the current species scope.

        Args:
            frame_name: One of ``pathway``, ``term2gene``, or ``term2name``.

        Returns:
            The requested lazy frame. Parsing is deferred until the first frame
            request, then the frame mapping is reused by subsequent operations.

        Raises:
            KeyError: If ``frame_name`` is not a declared tidy frame.
            ValueError: If the GMT content is malformed when first parsed.

        Examples:
            Materialize only the pathway frame:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.lazy_frame("pathway").collect().height
            2
        """
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

    Examples:
        Create a selection through a local dataset handle:

        >>> db = WikiPathwaysDb.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
        ... )
        >>> selection = db.select_ids(["2687"])
        >>> isinstance(selection, WikiPathwaysSelection)
        True
    """

    dataset: WikiPathwaysDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _lf_mapping: pl.LazyFrame | None = field(default=None, repr=False)
    _lf_unmapped: pl.LazyFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> selection = db.select_groups({"case": ["2687"]})
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the group ID column when this selection is grouped."""
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every pathway mapping matched by the selected Entrez IDs.

        Examples:
            Materialize the pathway mapped to Entrez ID 2687:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> selection = db.select_ids(["2687"])
            >>> selection.extract_mapping()["WikiPathwaysId"].to_list()
            ['WP100']
        """
        return self._lazy_mapping().collect()

    def _lazy_mapping(self) -> pl.LazyFrame:
        if self._lf_mapping is None:
            self._lf_mapping = extract_mapping_frame(
                self.dataset.lazy_frame("pathway"),
                self.dataset.lazy_frame("term2gene"),
                self._df_input_ids,
                cols_group_id=self._col_group_id,
            )
        return self._lf_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no WikiPathways mapping.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain an Entrez ID absent from the snapshot:

            >>> db = WikiPathwaysDb.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> selection = db.select_ids(["2687", "MISSING"])
            >>> selection.extract_unmapped_input_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
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
