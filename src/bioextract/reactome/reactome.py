from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
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
    extract_unmatched_ids_frame,
    filter_relation_frame,
    filter_species_frame,
    read_mapping_frame,
    read_pathway_frame,
    read_relation_frame,
)

__all__ = [
    "ReactomeDatabase",
]


@dataclass(frozen=True, slots=True)
class _ReactomeSnapshot:
    file_uniprot2reactome: Path | None = None
    file_pathways: Path | None = None
    file_relations: Path | None = None


@dataclass(slots=True)
class ReactomeDatabase:
    """Path-first access to local Reactome mapping snapshots.

    `ReactomeDatabase` is the public entrypoint for extracting Reactome annotation
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

        >>> db = ReactomeDatabase.from_files(
        ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
        ...     pathways="data/reactome/ReactomePathways.txt",
        ... )
        >>> (
        ...     db.with_species("Homo sapiens")
        ...     .select_ids(["P04637"])
        ...     .extract_mapping()["ReactomePathwayId"]
        ...     .to_list()
        ... )
        ['R-HSA-6798695', 'R-HSA-69563']
    """

    snapshot: _ReactomeSnapshot
    species: str | None = None
    _df_mapping_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_pathways_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_relations_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_pathways: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_relations: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_term2gene: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_term2name: pl.DataFrame | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_files(
        cls,
        *,
        uniprot_mapping: os.PathLike[str] | str | None = None,
        pathways: os.PathLike[str] | str | None = None,
        relations: os.PathLike[str] | str | None = None,
    ) -> ReactomeDatabase:
        """Create a dataset handle from local Reactome files.

        Args:
            uniprot_mapping: Path to `UniProt2Reactome.txt`.
            pathways: Path to `ReactomePathways.txt`.
            relations: Path to `ReactomePathwaysRelation.txt`.

        Returns:
            A dataset handle that can build whole-resource frames or selections.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If no files are provided.

        Examples:
            Open compact mapping and pathway fixtures:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
            ...     pathways="data/reactome/ReactomePathways.txt",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['mapping', 'pathway', 'term2gene', 'term2name']
        """
        if uniprot_mapping is None and pathways is None and relations is None:
            raise ValueError("At least one Reactome input file must be provided")
        file_uniprot2reactome = uniprot_mapping
        file_pathways = pathways
        file_relations = relations
        if file_uniprot2reactome is not None:
            file_uniprot2reactome = _validate_reactome_file(
                file_uniprot2reactome,
                label="Reactome UniProt2Reactome file",
            )
        if file_pathways is not None:
            file_pathways = _validate_reactome_file(
                file_pathways,
                label="Reactome pathways file",
            )
        if file_relations is not None:
            file_relations = _validate_reactome_file(
                file_relations,
                label="Reactome pathway relations file",
            )
        return cls(
            snapshot=_ReactomeSnapshot(
                file_uniprot2reactome=file_uniprot2reactome,
                file_pathways=file_pathways,
                file_relations=file_relations,
            ),
        )

    def with_species(self, species: str) -> ReactomeDatabase:
        """Create a species-scoped view of this Reactome snapshot.

        Args:
            species: Reactome species display name, matched exactly after
                trimming whitespace.

        Returns:
            A new dataset handle sharing the same file paths.

        Raises:
            ValueError: If the normalized species string is empty.

        Examples:
            Exclude pathways from other species:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> (
            ...     db.with_species("Homo sapiens")
            ...     .build_tidy()
            ...     .frames["mapping"]
            ...     .collect()["Species"]
            ...     .unique()
            ...     .to_list()
            ... )
            ['Homo sapiens']
        """
        species_normalized = str(species).strip()
        if not species_normalized:
            raise ValueError("Reactome species must be non-empty after normalization")
        return ReactomeDatabase(
            snapshot=self.snapshot,
            species=species_normalized,
        )

    def select_ids(self, ids: Iterable[str]) -> ReactomeSelection:
        """Create a single-query selection from UniProt accessions.

        Args:
            ids: Input UniProt accessions. Pipe-style UniProt values are
                normalized by the shared input normalizer.

        Returns:
            A selection that can extract pathway mappings and unmapped IDs.

        Examples:
            Normalize a UniProt pipe ID and retain an unmapped accession:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_ids(
            ...     ["sp|P04637|P53_HUMAN", "MISSING"]
            ... )
            >>> selection.extract_mapping()["InputId"].unique().to_list()
            ['P04637']
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> ReactomeSelection:
        """Create a grouped selection from multiple UniProt accession sets.

        Args:
            ids_by_group: Mapping from group label to input UniProt accessions.

        Returns:
            A grouped selection that carries `GroupId` through outputs.

        Raises:
            ValueError: If group IDs are invalid after normalization.

        Examples:
            Keep mapped and unmapped accessions in their original groups:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_groups(
            ...     {"tumor": ["P04637"], "control": ["MISSING"]}
            ... )
            >>> (
            ...     selection.extract_mapping()
            ...     .select("GroupId", "InputId")
            ...     .unique()
            ...     .to_dicts()
            ... )
            [{'GroupId': 'tumor', 'InputId': 'P04637'}]
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'GroupId': 'control', 'InputId': 'MISSING'}]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            _df_group_membership=grp_in_frames.df_group_membership,
        )

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract distinct Reactome-pathway-to-UniProt enrichment pairs.

        The current species scope is applied before pairs are deduplicated and
        sorted by pathway and accession.

        Raises:
            ValueError: If the UniProt-to-Reactome mapping file was not supplied.

        Examples:
            Extract human pathway-to-UniProt enrichment pairs:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... ).with_species("Homo sapiens")
            >>> db.extract_term2gene().head(2).to_dicts()
            [{'ReactomePathwayId': 'R-HSA-6798695', 'UniProtId': 'P04637'}, {'ReactomePathwayId': 'R-HSA-6798695', 'UniProtId': 'Q9Y243'}]
        """
        if self._df_term2gene is None:
            self._df_term2gene = extract_term2gene_frame(self._mapping_frame())
        return self._df_term2gene

    def extract_term2name(self) -> pl.DataFrame:
        """Extract pathway display names and species metadata for enrichment.

        Returns:
            One row per Reactome pathway ID in the current species scope.

        Raises:
            ValueError: If the pathway metadata file was not supplied.

        Examples:
            Extract human pathway labels for enrichment output:

            >>> db = ReactomeDatabase.from_files(
            ...     pathways="data/reactome/ReactomePathways.txt"
            ... ).with_species("Homo sapiens")
            >>> (
            ...     db.extract_term2name()
            ...     .select("ReactomePathwayId", "PathwayName")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'ReactomePathwayId': 'R-HSA-1640170', 'PathwayName': 'Cell Cycle'}]
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

            >>> db = ReactomeDatabase.from_files(
            ...     pathways="data/reactome/ReactomePathways.txt",
            ...     relations="data/reactome/ReactomePathwaysRelation.txt",
            ... ).with_species("Homo sapiens")
            >>> db.extract_pathway_relations().head(1).to_dicts()
            [{'ParentReactomePathwayId': 'R-HSA-1640170', 'ChildReactomePathwayId': 'R-HSA-6798695'}]
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

    def build_tidy(self) -> TidyDataset:
        """Build a lazy Reactome tidy dataset from the available source files.

        Returns:
            A `TidyDataset` containing only assets derivable from the provided
            raw files. Species scoping is reflected in every applicable frame.

        Examples:
            Build all five tidy frames when every resource is available:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
            ...     pathways="data/reactome/ReactomePathways.txt",
            ...     relations="data/reactome/ReactomePathwaysRelation.txt",
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

        return TidyDataset(
            frames={frame_name: frame.lazy() for frame_name, frame in frames.items()},
            source=self._tidy_sources(),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="reactome-mapping-files-v1",
            build_id_prefix="reactome-mapping",
            assets=tuple(assets),
            resource_name="reactome",
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish available Reactome relations as one DuckDB.

        Enrichment projections are omitted because they duplicate the canonical
        protein-pathway and pathway relations.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(Path(dir_out) / "reactome.duckdb")
            ...     result.tables
            ('protein_pathway',)
        """
        dataset = self.build_tidy()
        assets = tuple(
            asset
            for asset in dataset.assets
            if asset.frame_name not in {"term2gene", "term2name"}
        )
        canonical = TidyDataset(
            frames=dataset.frames,
            source=dataset.source,
            resource_schema_version=dataset.resource_schema_version,
            source_schema_profile=dataset.source_schema_profile,
            build_id_prefix=dataset.build_id_prefix,
            assets=assets,
            resource_name=dataset.resource_name,
            release_version=dataset.release_version,
        )
        return canonical.write_duckdb(
            path,
            table_names={
                "mapping": "protein_pathway",
                "relation": "pathway_relation",
            },
            if_exists=if_exists,
        )

    def _mapping_frame(self) -> pl.DataFrame:
        """Return and cache the normalized mapping in the current scope.

        Raises:
            ValueError: If the UniProt-to-Reactome mapping file was not supplied.

        Notes:
            Public extractors own the stable output contracts. Keep this raw
            normalized frame private so parser columns cannot become an
            accidental API.
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
                return self._mapping_frame()
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
                    logical_name="uniprot_mapping",
                    path=self.snapshot.file_uniprot2reactome,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        if self.snapshot.file_pathways is not None:
            sources.append(
                TidySource(
                    logical_name="pathways",
                    path=self.snapshot.file_pathways,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        if self.snapshot.file_relations is not None:
            sources.append(
                TidySource(
                    logical_name="relations",
                    path=self.snapshot.file_relations,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        return tuple(sources)


@dataclass(slots=True)
class ReactomeSelection:
    """Selection handle for single and grouped Reactome queries.

    Selections are created by :meth:`ReactomeDatabase.select_ids` or
    :meth:`ReactomeDatabase.select_groups`. Single selections return tables keyed by
    `InputId`; grouped selections prepend `GroupId`.

    Examples:
        Use a returned selection to materialize matched pathways:

        >>> db = ReactomeDatabase.from_files(
        ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
        ... )
        >>> selection = db.select_ids(["P04637"])
        >>> (
        ...     selection.extract_mapping()
        ...     .select("InputId", "ReactomePathwayId")
        ...     .head(1)
        ...     .to_dicts()
        ... )
        [{'InputId': 'P04637', 'ReactomePathwayId': 'R-HSA-6798695'}]
    """

    dataset: ReactomeDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_groups({"tumor": ["P04637"]})
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every pathway mapping matched by the selected UniProt IDs.

        Examples:
            Materialize the two fixture pathways mapped to TP53:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... ).with_species("Homo sapiens")
            >>> selection = db.select_ids(["P04637"])
            >>> selection.extract_mapping()["ReactomePathwayId"].to_list()
            ['R-HSA-6798695', 'R-HSA-69563']
        """
        if self._df_mapping is None:
            self._df_mapping = extract_mapping_frame(
                self.dataset._mapping_frame(),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                self._df_input_ids,
                df_group_membership=self._df_group_membership,
            )
        return self._df_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no Reactome pathway mapping.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain a normalized accession with no Reactome mapping:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_ids(["P04637", "MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                df_group_membership=self._df_group_membership,
            )
        return self._df_unmapped


def _validate_reactome_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    return file_path
