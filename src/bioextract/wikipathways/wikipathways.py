from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._publication import DuckDBWriteResult
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource

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
    extract_unmatched_ids_frame,
    read_gmt_frames,
    resolve_gmt_sources,
)

__all__ = [
    "WikiPathwaysDatabase",
    "WikiPathwaysTidyDataset",
]


@dataclass(frozen=True, slots=True)
class _WikiPathwaysSnapshot:
    files_gmt: tuple[Path, ...]
    species: str | None = None


WikiPathwaysTidyDataset = TidyDataset


@dataclass(slots=True)
class WikiPathwaysDatabase:
    """Source-first access to a local WikiPathways GMT snapshot.

    `WikiPathwaysDatabase` is the public entrypoint for extracting pathway gene sets
    from one or more local WikiPathways GMT files. GMT rows are interpreted as
    pathway metadata followed by NCBI Entrez Gene IDs; the class does not
    perform identifier conversion or enrichment statistics.

    Construct instances with :meth:`from_gmt`, then either extract whole
    resource frames or create single/grouped Entrez ID selections through
    :meth:`select_ids` and :meth:`select_groups`.

    Examples:
        Extract pathway mappings while preserving unmapped Entrez IDs:

        >>> db = WikiPathwaysDatabase.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
        ...     species="Homo sapiens",
        ... )
        >>> selection = db.select_ids(["2687", "MISSING"])
        >>> (
        ...     selection.extract_mapping()
        ...     .select("InputId", "WikiPathwaysId")
        ...     .to_dicts()
        ... )
        [{'InputId': '2687', 'WikiPathwaysId': 'WP100'}]
        >>> selection.extract_unmatched_ids().to_dicts()
        [{'InputId': 'MISSING'}]
    """

    snapshot: _WikiPathwaysSnapshot
    _frames: dict[str, pl.LazyFrame] | None = field(
        default=None, init=False, repr=False
    )
    _release_version: str | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_gmt(
        cls,
        source: os.PathLike[str] | str | Sequence[os.PathLike[str] | str],
        *,
        species: str | None = None,
        glob: bool = True,
    ) -> WikiPathwaysDatabase:
        """Create a dataset handle from one or more local WikiPathways GMT files.

        Args:
            source: One local path, a sequence of local paths, or glob
                expressions. Every resolved physical file must be unique.
            species: Optional species display name used as an exact metadata
                filter after parsing.
            glob: Expand each source entry as a glob expression. If false,
                treat every entry literally.

        Returns:
            A dataset handle that can build pathway, term2gene, and term2name
            frames.

        Raises:
            FileNotFoundError: If a literal file does not exist or a glob
                expression has no matches.
            ValueError: At construction, if the source is empty, resolves to a
                directory, repeats a physical file, or the species string is
                empty after normalization. GMT content parsing is deferred;
                malformed content, inconsistent Collection or Version values,
                and duplicate WikiPathways IDs raise when a frame is first
                extracted, built, or written.

        Examples:
            Open all GMT files in one local snapshot and retain human rows:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/*.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> (
            ...     db.extract_pathway()
            ...     .select("WikiPathwaysId", "PathwayName")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'WikiPathwaysId': 'WP100', 'PathwayName': 'Glutathione metabolism'}]
        """
        files_gmt = resolve_gmt_sources(source, glob=glob)
        species_normalized = None if species is None else str(species).strip()
        if species is not None and not species_normalized:
            raise ValueError(
                "WikiPathways species must be non-empty after normalization"
            )
        return cls(
            snapshot=_WikiPathwaysSnapshot(
                files_gmt=files_gmt,
                species=species_normalized,
            ),
        )

    def select_ids(self, ids: Iterable[str]) -> WikiPathwaysSelection:
        """Create a single-query selection from Entrez Gene IDs.

        Args:
            ids: Input Entrez Gene IDs. Values are normalized with the shared
                input-ID normalizer.

        Returns:
            A selection that can extract pathway mappings and unmapped IDs.

        Examples:
            Select the pathway containing Entrez Gene ID 2687:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.select_ids(["2687"])
            ...     .extract_mapping()
            ...     .select("InputId", "WikiPathwaysId")
            ...     .to_dicts()
            ... )
            [{'InputId': '2687', 'WikiPathwaysId': 'WP100'}]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return WikiPathwaysSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> WikiPathwaysSelection:
        """Create a grouped selection from multiple Entrez Gene ID sets.

        Args:
            ids_by_group: Mapping from group label to input Entrez Gene IDs.

        Returns:
            A grouped selection that carries `GroupId` through outputs.

        Raises:
            ValueError: If group IDs are invalid after normalization.

        Examples:
            Preserve comparison labels in the mapping table:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.select_groups({"case": ["2687"], "control": ["435"]})
            ...     .extract_mapping()
            ...     .select("GroupId", "InputId", "WikiPathwaysId")
            ...     .to_dicts()
            ... )
            [{'GroupId': 'case', 'InputId': '2687', 'WikiPathwaysId': 'WP100'}, {'GroupId': 'control', 'InputId': '435', 'WikiPathwaysId': 'WP106'}]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
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

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_pathway()["WikiPathwaysId"].to_list()
            ['WP100', 'WP106']
        """
        return self._lazy_frame("pathway").collect()

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract distinct WikiPathways-pathway-to-Entrez enrichment pairs.

        Examples:
            Inspect pathway-to-gene pairs from the compact human fixture:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_term2gene().head(2).to_dicts()
            [{'WikiPathwaysId': 'WP100', 'GeneId': '2678'}, {'WikiPathwaysId': 'WP100', 'GeneId': '2687'}]
        """
        return self._lazy_frame("term2gene").collect()

    def extract_term2name(self) -> pl.DataFrame:
        """Extract one pathway display-metadata row per WikiPathways ID.

        Examples:
            Look up a display name for an enrichment term:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.extract_term2name()
            ...     .select("WikiPathwaysId", "PathwayName")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'WikiPathwaysId': 'WP100', 'PathwayName': 'Glutathione metabolism'}]
        """
        return self._lazy_frame("term2name").collect()

    def build_tidy(self) -> WikiPathwaysTidyDataset:
        """Build the lazy WikiPathways tidy dataset for the current scope.

        Returns:
            A `TidyDataset` with species-consistent ``pathway``, ``term2gene``,
            and ``term2name`` frames.

        Examples:
            Build the three declared WikiPathways frames:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> sorted(db.build_tidy().frames)
            ['pathway', 'term2gene', 'term2name']
        """
        frames = {
            "pathway": self._lazy_frame("pathway"),
            "term2gene": self._lazy_frame("term2gene"),
            "term2name": self._lazy_frame("term2name"),
        }
        release_version = self._release_version
        if release_version is None:
            raise RuntimeError("WikiPathways GMT release Version was not parsed")
        return WikiPathwaysTidyDataset(
            frames=frames,
            source=tuple(
                TidySource(
                    logical_name=f"pathway_gmt_{index:03d}",
                    path=file_gmt,
                    media_type=MEDIA_TYPE_GMT,
                )
                for index, file_gmt in enumerate(
                    self.snapshot.files_gmt,
                    start=1,
                )
            ),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="wikipathways-gmt-v1",
            build_id_prefix="wikipathways-gmt",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
            resource_name="wikipathways",
            release_version=release_version,
            release_version_source="official_metadata",
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish pathway and membership relations as one DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(
            ...         Path(dir_out) / "wikipathways.duckdb"
            ...     )
            ...     result.tables
            ('pathway', 'pathway_gene')
        """
        dataset = self.build_tidy()
        canonical = WikiPathwaysTidyDataset(
            frames=dataset.frames,
            source=dataset.source,
            resource_schema_version=dataset.resource_schema_version,
            source_schema_profile=dataset.source_schema_profile,
            build_id_prefix=dataset.build_id_prefix,
            assets=tuple(
                asset for asset in dataset.assets if asset.frame_name != "term2name"
            ),
            resource_name=dataset.resource_name,
            release_version=dataset.release_version,
            release_version_source=dataset.release_version_source,
        )
        return canonical.write_duckdb(
            path,
            table_names={"term2gene": "pathway_gene"},
            if_exists=if_exists,
        )

    def _lazy_frame(self, frame_name: str) -> pl.LazyFrame:
        """Return a cached lazy tidy frame in the current species scope.

        Args:
            frame_name: One of ``pathway``, ``term2gene``, or ``term2name``.

        Returns:
            The requested lazy frame. Parsing is deferred until the first frame
            request, then the frame mapping is reused by subsequent operations.

        Raises:
            KeyError: If ``frame_name`` is not a declared tidy frame.
            ValueError: If the GMT content is malformed when first parsed.

        Notes:
            Public extractors and `build_tidy()` own the supported frame
            contracts. Keep this cache entrypoint private so internal frame
            keys do not become a second API.
        """
        if self._frames is None:
            parsed = read_gmt_frames(self.snapshot.files_gmt)
            frames = parsed.frames
            self._release_version = parsed.release_version
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

    Selections are created by :meth:`WikiPathwaysDatabase.select_ids` or
    :meth:`WikiPathwaysDatabase.select_groups`. Single selections return tables keyed
    by `InputId`; grouped selections prepend `GroupId`.

    Examples:
        Use a returned selection to materialize matched pathways:

        >>> db = WikiPathwaysDatabase.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
        ... )
        >>> selection = db.select_ids(["2687"])
        >>> (
        ...     selection.extract_mapping()
        ...     .select("InputId", "WikiPathwaysId")
        ...     .to_dicts()
        ... )
        [{'InputId': '2687', 'WikiPathwaysId': 'WP100'}]
    """

    dataset: WikiPathwaysDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _lf_mapping: pl.LazyFrame | None = field(default=None, repr=False)
    _lf_unmapped: pl.LazyFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = WikiPathwaysDatabase.from_gmt(
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

            >>> db = WikiPathwaysDatabase.from_gmt(
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
                self.dataset._lazy_frame("pathway"),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                self.dataset._lazy_frame("term2gene"),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                self._df_input_ids,
                cols_group_id=self._col_group_id,
            )
        return self._lf_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no WikiPathways mapping.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain an Entrez ID absent from the snapshot:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> selection = db.select_ids(["2687", "MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._lf_unmapped is None:
            self._lf_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self._lazy_mapping(),
                cols_group_id=self._col_group_id,
            )
        return self._lf_unmapped.collect()


def _filter_species_frame(lf: pl.LazyFrame, species: str) -> pl.LazyFrame:
    if "Species" not in lf.collect_schema().names():
        return lf
    return lf.filter(pl.col("Species") == species)
