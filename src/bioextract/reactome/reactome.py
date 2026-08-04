from __future__ import annotations

import os
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import polars as pl

from bioextract._publication import (
    BIOEXTRACT_RELATIONS,
    DuckDBWriteResult,
    validate_duckdb_metadata_v1,
)
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError, IntegrityError

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


class _ReopenedReactomeTidyDataset(TidyDataset):
    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        table_names: Mapping[str, str] | None = None,
        if_exists: str = "fail",
        preserve_source_headers: Collection[str] = (),
        include_source_hashes: bool = False,
    ) -> DuckDBWriteResult:
        del path, table_names, if_exists, preserve_source_headers, include_source_hashes
        raise CapabilityError("write_duckdb() requires a Reactome source-file handle")


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
    _publication_path: Path | None = field(default=None, init=False, repr=False)
    _publication_identity: tuple[int, int, int, int, int] | None = field(
        default=None, init=False, repr=False
    )
    _publication_tables: frozenset[str] = field(
        default=frozenset(), init=False, repr=False
    )

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

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> ReactomeDatabase:
        """Open a validated Reactome publication for domain and SQL access.

        Validation reads metadata and catalog schemas only; it does not recount
        the biological relations. The returned handle is pinned to the exact
        file that passed validation.

        Args:
            path: A bioextract Reactome metadata-v1 DuckDB publication.

        Returns:
            A publication-backed handle with the capabilities recorded by the
            publication.

        Raises:
            FileNotFoundError: If the publication does not exist.
            IntegrityError: If its metadata, capability inventory, or physical
                schema is incompatible, or if the file changes during validation.

        Examples:
            Reopen a publication and select one UniProt accession:

            >>> db = ReactomeDatabase.from_duckdb(  # doctest: +SKIP
            ...     "tidy/reactome.duckdb"
            ... )
            >>> db.select_ids(["P04637"]).extract_mapping().height > 0  # doctest: +SKIP
            True
        """
        publication_path = Path(path).absolute()
        identity_before = _file_identity(publication_path)
        try:
            tables = _validate_reactome_publication(publication_path)
            identity_after = _file_identity(publication_path)
        except (KeyError, TypeError, ValueError) as error:
            raise IntegrityError(str(error)) from error
        except OSError as error:
            raise IntegrityError(
                "Reactome publication changed during validation"
            ) from error
        if identity_after != identity_before:
            raise IntegrityError("Reactome publication changed during validation")
        result = cls(snapshot=_ReactomeSnapshot())
        result._publication_path = publication_path
        result._publication_identity = identity_after
        result._publication_tables = tables
        return result

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a fresh caller-owned read-only DuckDB connection.

        Raises:
            CapabilityError: If this handle was created from source files.
            IntegrityError: If the validated publication was replaced or became
                unavailable.

        Examples:
            Run native SQL against a reopened publication:

            >>> db = ReactomeDatabase.from_duckdb(  # doctest: +SKIP
            ...     "tidy/reactome.duckdb"
            ... )
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
            >>> count >= 0  # doctest: +SKIP
            True
        """
        path = self._publication_path
        if path is None:
            raise CapabilityError("connect() requires ReactomeDatabase.from_duckdb()")
        self._assert_publication_identity()
        try:
            connection = duckdb.connect(str(path), read_only=True)
        except duckdb.Error as error:
            raise IntegrityError(
                "Reactome publication became unavailable; reopen it with from_duckdb()"
            ) from error
        try:
            self._assert_publication_identity()
        except BaseException:
            connection.close()
            raise
        return connection

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
        result = ReactomeDatabase(
            snapshot=self.snapshot,
            species=species_normalized,
        )
        result._publication_path = self._publication_path
        result._publication_identity = self._publication_identity
        result._publication_tables = self._publication_tables
        return result

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
            >>> selection.extract_mapping()["input_id"].unique().to_list()
            ['P04637']
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'input_id': 'MISSING'}]
        """
        self._assert_publication_current()
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
            A grouped selection that carries `group_id` through outputs.

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
            ...     .select("group_id", "input_id")
            ...     .unique()
            ...     .to_dicts()
            ... )
            [{'group_id': 'tumor', 'input_id': 'P04637'}]
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'group_id': 'control', 'input_id': 'MISSING'}]
        """
        self._assert_publication_current()
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
        self._assert_publication_current()
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
        self._assert_publication_current()
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
        self._assert_publication_current()
        if not self._has_relation():
            self._raise_missing_capability(
                "Cannot extract Reactome relations without relations file",
                "Reactome publication does not contain pathway relations",
            )
        if self._df_relations is None:
            if not self._has_pathway():
                if self.species is not None:
                    self._raise_missing_capability(
                        "Cannot apply Reactome species-scoped relation filtering "
                        "without pathways file",
                        "Reactome publication cannot apply species-scoped relation "
                        "filtering without pathway metadata",
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
        self._assert_publication_current()
        frames: dict[str, pl.DataFrame] = {}
        assets: list[TidyAsset] = []
        for path, kind, frame_name in ASSET_SPECS:
            if frame_name == "mapping" and not self._has_mapping():
                continue
            if frame_name == "term2gene" and not self._has_mapping():
                continue
            if frame_name == "pathway" and not self._has_pathway():
                continue
            if frame_name == "term2name" and not self._has_pathway():
                continue
            if frame_name == "relation" and not self._has_relation():
                continue
            frames[frame_name] = self._build_tidy_frame(frame_name)
            assets.append(TidyAsset(path=path, kind=kind, frame_name=frame_name))

        dataset_type = (
            _ReopenedReactomeTidyDataset
            if self._publication_path is not None
            else TidyDataset
        )
        return dataset_type(
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
        if self._publication_path is not None:
            raise CapabilityError(
                "write_duckdb() requires a Reactome source-file handle"
            )
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
        if not self._has_mapping():
            self._raise_missing_capability(
                "Cannot extract Reactome mapping without UniProt2Reactome file",
                "Reactome publication does not contain protein-pathway mappings",
            )
        if self._df_mapping is None:
            self._df_mapping = filter_species_frame(
                self._mapping_raw_frame(),
                self.species,
            )
        return self._df_mapping

    def _pathway_frame(self) -> pl.DataFrame:
        if not self._has_pathway():
            self._raise_missing_capability(
                "Cannot extract Reactome pathways without pathways file",
                "Reactome publication does not contain pathway metadata",
            )
        if self._df_pathways is None:
            self._df_pathways = filter_species_frame(
                self._pathway_raw_frame(),
                self.species,
            )
        return self._df_pathways

    def _mapping_raw_frame(self) -> pl.DataFrame:
        if not self._has_mapping():
            self._raise_missing_capability(
                "Cannot read Reactome mapping without UniProt2Reactome file",
                "Reactome publication does not contain protein-pathway mappings",
            )
        if self._df_mapping_raw is None:
            if self._publication_path is None:
                assert self.snapshot.file_uniprot2reactome is not None
                self._df_mapping_raw = read_mapping_frame(
                    self.snapshot.file_uniprot2reactome
                )
            else:
                self._df_mapping_raw = self._read_publication_table(
                    "protein_pathway",
                    {
                        "uniprot_id": "UniProtId",
                        "reactome_pathway_id": "ReactomePathwayId",
                        "reactome_url": "ReactomeUrl",
                        "pathway_name": "PathwayName",
                        "evidence_code": "EvidenceCode",
                        "species": "Species",
                    },
                )
        return self._df_mapping_raw

    def _pathway_raw_frame(self) -> pl.DataFrame:
        if not self._has_pathway():
            self._raise_missing_capability(
                "Cannot read Reactome pathways without pathways file",
                "Reactome publication does not contain pathway metadata",
            )
        if self._df_pathways_raw is None:
            if self._publication_path is None:
                assert self.snapshot.file_pathways is not None
                self._df_pathways_raw = read_pathway_frame(self.snapshot.file_pathways)
            else:
                self._df_pathways_raw = self._read_publication_table(
                    "pathway",
                    {
                        "reactome_pathway_id": "ReactomePathwayId",
                        "pathway_name": "PathwayName",
                        "species": "Species",
                    },
                )
        return self._df_pathways_raw

    def _relation_raw_frame(self) -> pl.DataFrame:
        if not self._has_relation():
            self._raise_missing_capability(
                "Cannot read Reactome relations without relations file",
                "Reactome publication does not contain pathway relations",
            )
        if self._df_relations_raw is None:
            if self._publication_path is None:
                assert self.snapshot.file_relations is not None
                self._df_relations_raw = read_relation_frame(
                    self.snapshot.file_relations
                )
            else:
                self._df_relations_raw = self._read_publication_table(
                    "pathway_relation",
                    {
                        "parent_reactome_pathway_id": "ParentReactomePathwayId",
                        "child_reactome_pathway_id": "ChildReactomePathwayId",
                    },
                )
        return self._df_relations_raw

    def _has_mapping(self) -> bool:
        return (
            self.snapshot.file_uniprot2reactome is not None
            or "protein_pathway" in self._publication_tables
        )

    def _has_pathway(self) -> bool:
        return (
            self.snapshot.file_pathways is not None
            or "pathway" in self._publication_tables
        )

    def _has_relation(self) -> bool:
        return (
            self.snapshot.file_relations is not None
            or "pathway_relation" in self._publication_tables
        )

    def _assert_publication_identity(self) -> None:
        path = self._publication_path
        try:
            current_identity = None if path is None else _file_identity(path)
        except OSError:
            current_identity = None
        if current_identity != self._publication_identity:
            raise IntegrityError(
                "Reactome publication was replaced; reopen it with from_duckdb()"
            )

    def _assert_publication_current(self) -> None:
        if self._publication_path is not None:
            self._assert_publication_identity()

    def _raise_missing_capability(
        self,
        source_message: str,
        publication_message: str,
    ) -> None:
        if self._publication_path is not None:
            raise CapabilityError(publication_message)
        raise ValueError(source_message)

    def _read_publication_table(
        self,
        table_name: str,
        rename: Mapping[str, str],
    ) -> pl.DataFrame:
        with self.connect() as connection:
            frame = pl.read_database(  # pyright: ignore[reportUnknownMemberType]
                f'SELECT * FROM "{table_name}"', connection
            )
        return frame.rename(rename)

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
    `input_id`; grouped selections prepend `group_id`.

    Examples:
        Use a returned selection to materialize matched pathways:

        >>> db = ReactomeDatabase.from_files(
        ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
        ... )
        >>> selection = db.select_ids(["P04637"])
        >>> (
        ...     selection.extract_mapping()
        ...     .select("input_id", "ReactomePathwayId")
        ...     .head(1)
        ...     .to_dicts()
        ... )
        [{'input_id': 'P04637', 'ReactomePathwayId': 'R-HSA-6798695'}]
    """

    dataset: ReactomeDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `group_id` through outputs.

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
        self.dataset._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
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
        group and include ``group_id`` in the result.

        Examples:
            Retain a normalized accession with no Reactome mapping:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> selection = db.select_ids(["P04637", "MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'input_id': 'MISSING'}]
        """
        self.dataset._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
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


_REACTOME_TABLE_CONTRACTS: dict[str, tuple[str, str, tuple[tuple[str, str], ...]]] = {
    "protein_pathway": (
        "uniprot_mapping",
        "canonical",
        (
            ("uniprot_id", "VARCHAR"),
            ("reactome_pathway_id", "VARCHAR"),
            ("reactome_url", "VARCHAR"),
            ("pathway_name", "VARCHAR"),
            ("evidence_code", "VARCHAR"),
            ("species", "VARCHAR"),
        ),
    ),
    "pathway": (
        "pathways",
        "canonical",
        (
            ("reactome_pathway_id", "VARCHAR"),
            ("pathway_name", "VARCHAR"),
            ("species", "VARCHAR"),
        ),
    ),
    "pathway_relation": (
        "relations",
        "canonical",
        (
            ("parent_reactome_pathway_id", "VARCHAR"),
            ("child_reactome_pathway_id", "VARCHAR"),
        ),
    ),
}

_REACTOME_COLUMN_MAPPINGS = {
    ("pathway", "PathwayName", "pathway_name", "generated_snake_case"),
    ("pathway", "ReactomePathwayId", "reactome_pathway_id", "generated_snake_case"),
    ("pathway", "Species", "species", "generated_snake_case"),
    (
        "pathway_relation",
        "ChildReactomePathwayId",
        "child_reactome_pathway_id",
        "generated_snake_case",
    ),
    (
        "pathway_relation",
        "ParentReactomePathwayId",
        "parent_reactome_pathway_id",
        "generated_snake_case",
    ),
    ("protein_pathway", "EvidenceCode", "evidence_code", "generated_snake_case"),
    ("protein_pathway", "PathwayName", "pathway_name", "generated_snake_case"),
    (
        "protein_pathway",
        "ReactomePathwayId",
        "reactome_pathway_id",
        "generated_snake_case",
    ),
    ("protein_pathway", "ReactomeUrl", "reactome_url", "generated_snake_case"),
    ("protein_pathway", "Species", "species", "generated_snake_case"),
    ("protein_pathway", "UniProtId", "uniprot_id", "generated_snake_case"),
}


def _validate_reactome_publication(path: Path) -> frozenset[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with duckdb.connect(str(path), read_only=True) as connection:
            metadata_rows = connection.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
            metadata = {str(row[0]): str(row[1]) for row in metadata_rows}
            if len(metadata) != len(metadata_rows):
                raise ValueError("Reactome publication has duplicate metadata keys")
            if metadata.get("bioextract.metadata_schema_version") != "1":
                raise ValueError("Unsupported Reactome metadata schema version")
            validate_duckdb_metadata_v1(connection, metadata)
            if metadata.get("bioextract.resource_name") != "reactome":
                raise ValueError("DuckDB file is not a bioextract Reactome publication")
            if metadata.get("bioextract.source_schema_profile") != (
                "reactome-mapping-files-v1"
            ):
                raise ValueError("Unsupported Reactome source schema profile")
            if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported Reactome resource schema version")
            if "bioextract.scope" in metadata:
                raise ValueError("Reactome publication scope is unsupported")

            source_rows = connection.execute(
                "SELECT logical_name, bytes FROM _bioextract.source_file"
            ).fetchall()
            source_roles = {str(row[0]) for row in source_rows}
            allowed_roles = {
                contract[0] for contract in _REACTOME_TABLE_CONTRACTS.values()
            }
            if (
                not source_roles
                or len(source_roles) != len(source_rows)
                or not source_roles <= allowed_roles
                or any(int(row[1]) < 0 for row in source_rows)
            ):
                raise ValueError("Reactome source capability inventory is unsupported")

            expected_tables = {
                table_name
                for table_name, (source_role, _role, _schema) in (
                    _REACTOME_TABLE_CONTRACTS.items()
                )
                if source_role in source_roles
            }
            relations = {
                (str(row[0]), str(row[1]), str(row[2]))
                for row in connection.execute(
                    "SELECT table_schema, table_name, table_type "
                    "FROM information_schema.tables "
                    "WHERE table_schema NOT IN ('information_schema', 'pg_catalog')"
                ).fetchall()
            }
            expected_relations = {
                ("_bioextract", name, "BASE TABLE") for name in BIOEXTRACT_RELATIONS
            } | {("main", name, "BASE TABLE") for name in expected_tables}
            if relations != expected_relations:
                raise ValueError(
                    "Reactome physical table/view inventory is unsupported"
                )

            info_rows = connection.execute(
                "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
            ).fetchall()
            recorded = {str(row[0]): (str(row[1]), int(row[2])) for row in info_rows}
            if len(recorded) != len(info_rows) or set(recorded) != expected_tables:
                raise ValueError("Reactome table inventory does not match metadata")
            for table_name, (role, row_count) in recorded.items():
                _source_role, expected_role, expected_schema = (
                    _REACTOME_TABLE_CONTRACTS[table_name]
                )
                if role != expected_role or row_count < 0:
                    raise ValueError(
                        "Reactome table capability inventory is unsupported"
                    )
                actual_schema = tuple(
                    (str(row[1]), str(row[2]))
                    for row in connection.execute(
                        f'PRAGMA table_info("{table_name}")'
                    ).fetchall()
                )
                if actual_schema != expected_schema:
                    raise ValueError(
                        f"Reactome table schema is unsupported: {table_name}"
                    )
            column_mappings = {
                tuple(str(value) for value in row)
                for row in connection.execute(
                    "SELECT table_name, source_column, output_column, reason "
                    "FROM _bioextract.column_mapping"
                ).fetchall()
            }
            expected_column_mappings = {
                row for row in _REACTOME_COLUMN_MAPPINGS if row[0] in expected_tables
            }
            if column_mappings != expected_column_mappings:
                raise ValueError("Reactome column provenance inventory is unsupported")
    except duckdb.Error as error:
        raise ValueError(f"Cannot open Reactome DuckDB publication: {path}") from error
    return frozenset(expected_tables)


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
