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
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    extract_mapping_frame,
    extract_term2gene_frame,
    extract_term2name_frame,
    extract_unmapped_input_ids_frame,
    filter_relation_frame,
    filter_species_frame,
    read_mapping_frame,
    read_pathway_frame,
    read_relation_frame,
)

__all__ = [
    "ReactomeDb",
    "ReactomeResourceLimits",
    "ReactomeTidyDataset",
]


@dataclass(frozen=True, slots=True)
class ReactomeResourceLimits:
    """Optional fail-fast limits for Reactome resources and selections.

    File limits are measured in bytes and checked when a snapshot handle is
    created. Count limits apply after input IDs and group names are normalized;
    ``None`` disables the corresponding limit.

    Examples:
        Limit grouped selections to 20 groups:

        >>> limits = ReactomeResourceLimits(num_groups_max=20)
        >>> limits.num_groups_max
        20
    """

    file_uniprot2reactome_bytes_max: int | None = None
    file_pathways_bytes_max: int | None = None
    file_relations_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _ReactomeSnapshot:
    file_uniprot2reactome: Path | None = None
    file_pathways: Path | None = None
    file_relations: Path | None = None


ReactomeTidyDataset = TidyDataset


@dataclass(slots=True)
class ReactomeDb:
    """Path-first access to local Reactome mapping snapshots.

    `ReactomeDb` is the public entrypoint for extracting Reactome annotation
    mappings and standard enrichment inputs from local open-data files. The
    three raw files are composable: callers may provide only the files needed
    by the requested capability, and missing-file errors are raised at the
    feature boundary.

    Construct instances with :meth:`from_files`, optionally constrain them with
    :meth:`with_species`, then either extract whole-resource frames or create
    single/grouped selections through :meth:`select_ids` and
    :meth:`select_groups`.

    Examples:
        Select UniProt mappings from one species:

        >>> db = ReactomeDb.from_files(
        ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt",
        ...     file_pathways="data/reactome/ReactomePathways.txt",
        ... )
        >>> db.with_species("Homo sapiens").select_ids(
        ...     ["P04637"]
        ... ).extract_mapping().height
        2
    """

    snapshot: _ReactomeSnapshot
    limits: ReactomeResourceLimits = field(default_factory=ReactomeResourceLimits)
    species: str | None = None
    _df_mapping_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_pathways_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_relations_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_pathways: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_relations: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_term2gene: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_term2name: pl.DataFrame | None = field(default=None, init=False, repr=False)

    DEFAULT_RESOURCE_LIMITS = ReactomeResourceLimits()

    @classmethod
    def from_files(
        cls,
        *,
        file_uniprot2reactome: os.PathLike[str] | str | None = None,
        file_pathways: os.PathLike[str] | str | None = None,
        file_relations: os.PathLike[str] | str | None = None,
        limits: ReactomeResourceLimits | None = None,
    ) -> ReactomeDb:
        """Create a dataset handle from local Reactome files.

        Args:
            file_uniprot2reactome: Path to `UniProt2Reactome.txt`.
            file_pathways: Path to `ReactomePathways.txt`.
            file_relations: Path to `ReactomePathwaysRelation.txt`.
            limits: Optional resource policy. When omitted, no finite size or
                selection limits are imposed.

        Returns:
            A dataset handle that can build whole-resource frames or selections.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If no files are provided or a configured file-size limit
                is exceeded.

        Examples:
            Open compact mapping and pathway fixtures:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt",
            ...     file_pathways="data/reactome/ReactomePathways.txt",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['mapping', 'pathway', 'term2gene', 'term2name']
        """
        limits_resolved = ReactomeResourceLimits() if limits is None else limits
        if (
            file_uniprot2reactome is None
            and file_pathways is None
            and file_relations is None
        ):
            raise ValueError("At least one Reactome input file must be provided")
        if file_uniprot2reactome is not None:
            file_uniprot2reactome = _validate_reactome_file(
                file_uniprot2reactome,
                size_max=limits_resolved.file_uniprot2reactome_bytes_max,
                label="Reactome UniProt2Reactome file",
            )
        if file_pathways is not None:
            file_pathways = _validate_reactome_file(
                file_pathways,
                size_max=limits_resolved.file_pathways_bytes_max,
                label="Reactome pathways file",
            )
        if file_relations is not None:
            file_relations = _validate_reactome_file(
                file_relations,
                size_max=limits_resolved.file_relations_bytes_max,
                label="Reactome pathway relations file",
            )
        return cls(
            snapshot=_ReactomeSnapshot(
                file_uniprot2reactome=file_uniprot2reactome,
                file_pathways=file_pathways,
                file_relations=file_relations,
            ),
            limits=limits_resolved,
        )

    def with_species(self, species: str) -> ReactomeDb:
        """Create a species-scoped view of this Reactome snapshot.

        Args:
            species: Reactome species display name, matched exactly after
                trimming whitespace.

        Returns:
            A new dataset handle sharing the same file paths and limits.

        Raises:
            ValueError: If the normalized species string is empty.

        Examples:
            Create a human-only view without mutating the original handle:

            >>> db = ReactomeDb.from_files(
            ...     file_pathways="data/reactome/ReactomePathways.txt"
            ... )
            >>> db.with_species("Homo sapiens").species
            'Homo sapiens'
        """
        species_normalized = str(species).strip()
        if not species_normalized:
            raise ValueError("Reactome species must be non-empty after normalization")
        return ReactomeDb(
            snapshot=self.snapshot,
            limits=self.limits,
            species=species_normalized,
        )

    def select_ids(self, ids: Iterable[str]) -> ReactomeSelection:
        """Create a single-query selection from UniProt accessions.

        Args:
            ids: Input UniProt accessions. Pipe-style UniProt values are
                normalized by the shared input normalizer.

        Returns:
            A selection that can extract pathway mappings and unmapped IDs.

        Raises:
            ValueError: If the normalized input-ID count exceeds the configured
                limit.

        Examples:
            Create an ungrouped UniProt selection:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> db.select_ids(["P04637"]).is_grouped
            False
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
    ) -> ReactomeSelection:
        """Create a grouped selection from multiple UniProt accession sets.

        Args:
            group_to_ids: Mapping from group label to input UniProt accessions.

        Returns:
            A grouped selection that carries `GroupId` through outputs.

        Raises:
            ValueError: If group or input-ID limits are exceeded, or if group
                IDs are invalid after normalization.

        Examples:
            Create a grouped UniProt selection:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> db.select_groups({"tumor": ["P04637"]}).is_grouped
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
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
        )

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract distinct Reactome-pathway-to-UniProt enrichment pairs.

        The current species scope is applied before pairs are deduplicated and
        sorted by pathway and accession.

        Raises:
            ValueError: If the UniProt-to-Reactome mapping file was not supplied.

        Examples:
            Extract human pathway-to-UniProt enrichment pairs:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... ).with_species("Homo sapiens")
            >>> db.extract_term2gene().shape
            (3, 2)
        """
        if self._df_term2gene is None:
            self._df_term2gene = extract_term2gene_frame(self.mapping_frame())
        return self._df_term2gene

    def extract_term2name(self) -> pl.DataFrame:
        """Extract pathway display names and species metadata for enrichment.

        Returns:
            One row per Reactome pathway ID in the current species scope.

        Raises:
            ValueError: If the pathway metadata file was not supplied.

        Examples:
            Extract human pathway labels for enrichment output:

            >>> db = ReactomeDb.from_files(
            ...     file_pathways="data/reactome/ReactomePathways.txt"
            ... ).with_species("Homo sapiens")
            >>> db.extract_term2name().shape
            (3, 3)
        """
        if self._df_term2name is None:
            self._df_term2name = extract_term2name_frame(self._pathway_frame())
        return self._df_term2name

    def extract_pathway_relations(self) -> pl.DataFrame:
        """Extract Reactome parent-child pathway relations.

        When the dataset is species-scoped, both endpoints must exist in the
        species-scoped pathway metadata, so the pathways file is required.

        Raises:
            ValueError: If the relations file is missing, or if species-scoped
                filtering is requested without pathway metadata.

        Examples:
            Extract relations whose endpoints are both human pathways:

            >>> db = ReactomeDb.from_files(
            ...     file_pathways="data/reactome/ReactomePathways.txt",
            ...     file_relations="data/reactome/ReactomePathwaysRelation.txt",
            ... ).with_species("Homo sapiens")
            >>> db.extract_pathway_relations().height
            2
        """
        if self.snapshot.file_relations is None:
            raise ValueError("Cannot extract Reactome relations without relations file")
        if self._df_relations is None:
            if self.snapshot.file_pathways is None:
                if self.species is not None:
                    raise ValueError(
                        "Cannot apply Reactome species-scoped relation filtering "
                        "without pathways file"
                    )
                self._df_relations = (
                    self._relation_raw_frame()
                    .unique()
                    .sort(
                        "ParentReactomePathwayId",
                        "ChildReactomePathwayId",
                    )
                )
            else:
                self._df_relations = filter_relation_frame(
                    self._relation_raw_frame(),
                    self._pathway_frame(),
                )
        return self._df_relations

    def build_tidy(self) -> ReactomeTidyDataset:
        """Build a lazy Reactome tidy dataset from the available source files.

        Returns:
            A `TidyDataset` containing only assets derivable from the provided
            raw files. Species scoping is reflected in every applicable frame.

        Examples:
            Build all five tidy frames when every resource is available:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt",
            ...     file_pathways="data/reactome/ReactomePathways.txt",
            ...     file_relations="data/reactome/ReactomePathwaysRelation.txt",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['mapping', 'pathway', 'relation', 'term2gene', 'term2name']
        """
        frames: dict[str, pl.DataFrame] = {}
        assets: list[TidyAsset] = []
        for path, kind, frame_name in ASSET_SPECS:
            if frame_name == "mapping" and self.snapshot.file_uniprot2reactome is None:
                continue
            if (
                frame_name == "term2gene"
                and self.snapshot.file_uniprot2reactome is None
            ):
                continue
            if frame_name == "pathway" and self.snapshot.file_pathways is None:
                continue
            if frame_name == "term2name" and self.snapshot.file_pathways is None:
                continue
            if frame_name == "relation" and self.snapshot.file_relations is None:
                continue
            frames[frame_name] = self._build_tidy_frame(frame_name)
            assets.append(TidyAsset(path=path, kind=kind, frame_name=frame_name))

        return ReactomeTidyDataset(
            frames={frame_name: frame.lazy() for frame_name, frame in frames.items()},
            source=self._tidy_sources(),
            schema_version=SCHEMA_VERSION,
            build_id_prefix="reactome-mapping",
            assets=tuple(assets),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the Reactome tidy dataset as flat parquet files.

        Args:
            dir_out: Output directory for parquet assets.
            should_write_manifest: Whether to write `manifest.json`.
            should_hash_assets: Whether to calculate asset checksums in the
                manifest. This has no effect unless a manifest is requested.

        Returns:
            A write report with asset paths and optional manifest content.

        Examples:
            Write every asset available from a complete snapshot:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt",
            ...     file_pathways="data/reactome/ReactomePathways.txt",
            ...     file_relations="data/reactome/ReactomePathwaysRelation.txt",
            ... )
            >>> report = db.write_tidy("build/reactome")
            >>> [asset.path for asset in report.assets]
            ['mapping.parquet', 'pathway.parquet', 'relation.parquet', 'term2gene.parquet', 'term2name.parquet']
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def mapping_frame(self) -> pl.DataFrame:
        """Extract the raw UniProt-to-Reactome mapping in the current scope.

        Raises:
            ValueError: If the UniProt-to-Reactome mapping file was not supplied.

        Examples:
            Materialize the species-scoped raw mapping table:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... ).with_species("Homo sapiens")
            >>> db.mapping_frame().shape
            (3, 6)
        """
        if self.snapshot.file_uniprot2reactome is None:
            raise ValueError(
                "Cannot extract Reactome mapping without UniProt2Reactome file"
            )
        if self._df_mapping is None:
            self._df_mapping = filter_species_frame(
                self._mapping_raw_frame(),
                self.species,
            )
        return self._df_mapping

    def _pathway_frame(self) -> pl.DataFrame:
        if self.snapshot.file_pathways is None:
            raise ValueError("Cannot extract Reactome pathways without pathways file")
        if self._df_pathways is None:
            self._df_pathways = filter_species_frame(
                self._pathway_raw_frame(),
                self.species,
            )
        return self._df_pathways

    def _mapping_raw_frame(self) -> pl.DataFrame:
        if self.snapshot.file_uniprot2reactome is None:
            raise ValueError(
                "Cannot read Reactome mapping without UniProt2Reactome file"
            )
        if self._df_mapping_raw is None:
            self._df_mapping_raw = read_mapping_frame(
                self.snapshot.file_uniprot2reactome
            )
        return self._df_mapping_raw

    def _pathway_raw_frame(self) -> pl.DataFrame:
        if self.snapshot.file_pathways is None:
            raise ValueError("Cannot read Reactome pathways without pathways file")
        if self._df_pathways_raw is None:
            self._df_pathways_raw = read_pathway_frame(self.snapshot.file_pathways)
        return self._df_pathways_raw

    def _relation_raw_frame(self) -> pl.DataFrame:
        if self.snapshot.file_relations is None:
            raise ValueError("Cannot read Reactome relations without relations file")
        if self._df_relations_raw is None:
            self._df_relations_raw = read_relation_frame(self.snapshot.file_relations)
        return self._df_relations_raw

    def _build_tidy_frame(self, frame_name: str) -> pl.DataFrame:
        match frame_name:
            case "mapping":
                return self.mapping_frame()
            case "pathway":
                return self._pathway_frame()
            case "relation":
                return self.extract_pathway_relations()
            case "term2gene":
                return self.extract_term2gene()
            case "term2name":
                return self.extract_term2name()
            case _:
                raise ValueError(f"Unsupported Reactome tidy frame: {frame_name}")

    def _tidy_sources(self) -> tuple[TidySource, ...]:
        sources: list[TidySource] = []
        if self.snapshot.file_uniprot2reactome is not None:
            sources.append(
                TidySource(
                    path=self.snapshot.file_uniprot2reactome,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        if self.snapshot.file_pathways is not None:
            sources.append(
                TidySource(path=self.snapshot.file_pathways, media_type=MEDIA_TYPE_TSV)
            )
        if self.snapshot.file_relations is not None:
            sources.append(
                TidySource(path=self.snapshot.file_relations, media_type=MEDIA_TYPE_TSV)
            )
        return tuple(sources)


@dataclass(slots=True)
class ReactomeSelection:
    """Selection handle for single and grouped Reactome queries.

    Selections are created by :meth:`ReactomeDb.select_ids` or
    :meth:`ReactomeDb.select_groups`. Single selections return tables keyed by
    `InputId`; grouped selections prepend `GroupId`.

    Examples:
        Create a selection through a local dataset handle:

        >>> db = ReactomeDb.from_files(
        ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
        ... )
        >>> selection = db.select_ids(["P04637"])
        >>> isinstance(selection, ReactomeSelection)
        True
    """

    dataset: ReactomeDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_groups({"tumor": ["P04637"]})
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the group ID column when this selection is grouped."""
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every pathway mapping matched by the selected UniProt IDs.

        Examples:
            Materialize the two fixture pathways mapped to TP53:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... ).with_species("Homo sapiens")
            >>> selection = db.select_ids(["P04637"])
            >>> selection.extract_mapping()["ReactomePathwayId"].to_list()
            ['R-HSA-6798695', 'R-HSA-69563']
        """
        if self._df_mapping is None:
            self._df_mapping = extract_mapping_frame(
                self.dataset.mapping_frame(),
                self._df_input_ids,
                cols_group_id=self._col_group_id,
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no Reactome pathway mapping.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain a normalized accession with no Reactome mapping:

            >>> db = ReactomeDb.from_files(
            ...     file_uniprot2reactome="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_ids(["P04637", "MISSING"])
            >>> selection.extract_unmapped_input_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmapped_input_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                cols_group_id=self._col_group_id,
            )
        return self._df_unmapped


def _validate_reactome_file(
    file_path: os.PathLike[str] | str,
    *,
    size_max: int | None,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    validate_file_size(file_path=file_path, size_max=size_max, label=label)
    return file_path
