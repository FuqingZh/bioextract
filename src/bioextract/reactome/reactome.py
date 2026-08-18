from __future__ import annotations

import copy
import os
import re
from collections.abc import Collection, Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import (
    _RelationScanRequest,  # pyright: ignore[reportPrivateUsage]  # typed source boundary
    register_replayable_source,
)
from bioextract._publication import (
    BIOEXTRACT_RELATIONS,
    METADATA_SCHEMA_VERSION,
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    ValidationIssue,
    validate_duckdb_metadata_v2,
    write_duckdb_publication,
)
from bioextract._shared import (
    GroupInputFrames,
    normalize_input_id,
    validate_group_ids,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError, IntegrityError

from .constant import (
    ASSET_SPECS,
    COMPLEX_PATHWAY_ROLE,
    ENTITY_COLUMN_MAPPING_REASON,
    ENTITY_ROLE_BY_ROLE,
    EWAS_PATHWAY_ROLE,
    GMT_SOURCE_SPEC,
    MAPPING_OFFICIAL_FILENAMES,
    MAPPING_ROLE_BY_ARGUMENT,
    MAPPING_ROLE_BY_DIMENSIONS,
    MAPPING_ROLE_BY_ROLE,
    MAPPING_ROLE_SPECS,
    MEDIA_TYPE_TSV,
    MEDIA_TYPE_ZIP,
    PATHWAY_GENE_SET_ROLE,
    PATHWAY_ROLE,
    RELATION_ROLE,
    SCHEMA_VERSION,
    SOURCE_SCHEMA_PROFILE,
)
from .util import (
    extract_term2gene_frame,
    extract_term2name_frame,
    filter_relation_frame,
    filter_species_frame,
    read_entity_pathway_frame,
    read_gmt_frame,
    read_mapping_role_frame,
    read_pathway_frame,
    read_relation_frame,
    scan_entity_pathway_frame,
    scan_mapping_role_frame,
    scan_pathway_frame,
    scan_relation_frame,
)

__all__ = [
    "ReactomeDatabase",
]


@dataclass(frozen=True, slots=True)
class _ReactomeSnapshot:
    file_uniprot2reactome: Path | None = None
    file_uniprot_all_levels: Path | None = None
    file_uniprot_reactions: Path | None = None
    file_ncbi_mapping: Path | None = None
    file_ncbi_all_levels: Path | None = None
    file_ncbi_reactions: Path | None = None
    file_chebi_mapping: Path | None = None
    file_chebi_all_levels: Path | None = None
    file_chebi_reactions: Path | None = None
    file_gtop_mapping: Path | None = None
    file_gtop_all_levels: Path | None = None
    file_gtop_reactions: Path | None = None
    file_complex_pathways: Path | None = None
    file_ewas_pathways: Path | None = None
    file_pathway_gene_sets: Path | None = None
    file_pathways: Path | None = None
    file_relations: Path | None = None


class _ReopenedReactomeTidyDataset(TidyDataset):
    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        table_names: Mapping[str, str] | None = None,
        if_exists: str = "fail",
        source_columns: Mapping[str, Collection[str]] | None = None,
    ) -> DuckDBWriteResult:
        del path, table_names, if_exists, source_columns
        raise CapabilityError("write_duckdb() requires a Reactome source-file handle")


def _empty_frame_map() -> dict[str, pl.DataFrame]:
    return {}


@dataclass(slots=True)
class ReactomeDatabase:
    """Path-first access to local Reactome mapping snapshots.

    `ReactomeDatabase` is the public entrypoint for extracting Reactome annotation
    mappings and standard enrichment inputs from local open-data files. The
    The source roles are composable: callers may provide only the files needed
    by the requested capability, and missing-file errors are raised at the
    feature boundary. Every mapping namespace, target, pathway level, entity
    relation, and GMT relation remains a separate capability throughout
    construction, query, and publication.

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
        ...     .mappings().select("reactome_pathway_id").collect()
        ...     .to_series().to_list()
        ... )
        ['R-HSA-6798695', 'R-HSA-69563']
    """

    snapshot: _ReactomeSnapshot
    species: str | None = None
    _release_version: str | None = field(default=None, init=False, repr=False)
    _release_version_source: str | None = field(default=None, init=False, repr=False)
    _df_mapping_raw_by_role: dict[str, pl.DataFrame] = field(
        default_factory=_empty_frame_map, init=False, repr=False
    )
    _df_mapping_by_role: dict[str, pl.DataFrame] = field(
        default_factory=_empty_frame_map, init=False, repr=False
    )
    _df_entity_by_role: dict[str, pl.DataFrame] = field(
        default_factory=_empty_frame_map, init=False, repr=False
    )
    _df_pathway_gene_sets: pl.DataFrame | None = field(
        default=None, init=False, repr=False
    )
    _df_mapping_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_mapping_all_levels_raw: pl.DataFrame | None = field(
        default=None, init=False, repr=False
    )
    _df_pathways_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_relations_raw: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _df_mapping_all_levels: pl.DataFrame | None = field(
        default=None, init=False, repr=False
    )
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
        uniprot_all_levels: os.PathLike[str] | str | None = None,
        uniprot_reactions: os.PathLike[str] | str | None = None,
        ncbi_mapping: os.PathLike[str] | str | None = None,
        ncbi_all_levels: os.PathLike[str] | str | None = None,
        ncbi_reactions: os.PathLike[str] | str | None = None,
        chebi_mapping: os.PathLike[str] | str | None = None,
        chebi_all_levels: os.PathLike[str] | str | None = None,
        chebi_reactions: os.PathLike[str] | str | None = None,
        gtop_mapping: os.PathLike[str] | str | None = None,
        gtop_all_levels: os.PathLike[str] | str | None = None,
        gtop_reactions: os.PathLike[str] | str | None = None,
        complex_pathways: os.PathLike[str] | str | None = None,
        ewas_pathways: os.PathLike[str] | str | None = None,
        pathway_gene_sets: os.PathLike[str] | str | None = None,
        pathways: os.PathLike[str] | str | None = None,
        relations: os.PathLike[str] | str | None = None,
        release_version: str | None = None,
    ) -> ReactomeDatabase:
        """Create a dataset handle from local Reactome files.

        Args:
            uniprot_mapping: Path to the lowest-level
                `UniProt2Reactome.txt` relation.
            uniprot_all_levels: Path to the explicit hierarchy-expanded
                `UniProt2Reactome_All_Levels.txt` relation.
            pathways: Path to `ReactomePathways.txt`.
            relations: Path to `ReactomePathwaysRelation.txt`.
            release_version: Optional caller-declared Reactome snapshot
                identity. It is never inferred from a path or manifest.

        Returns:
            A dataset handle that can build whole-resource frames or selections.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If no files are provided.

        Examples:
            Open compact mapping and pathway fixtures:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
            ...     uniprot_all_levels="data/reactome/UniProt2Reactome_All_Levels.txt",
            ...     pathways="data/reactome/ReactomePathways.txt",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['pathway', 'uniprot_pathway_all_level', 'uniprot_pathway_lowest_level']
        """
        mapping_inputs = {
            name: value
            for name, value in {
                "uniprot_mapping": uniprot_mapping,
                "uniprot_all_levels": uniprot_all_levels,
                "uniprot_reactions": uniprot_reactions,
                "ncbi_mapping": ncbi_mapping,
                "ncbi_all_levels": ncbi_all_levels,
                "ncbi_reactions": ncbi_reactions,
                "chebi_mapping": chebi_mapping,
                "chebi_all_levels": chebi_all_levels,
                "chebi_reactions": chebi_reactions,
                "gtop_mapping": gtop_mapping,
                "gtop_all_levels": gtop_all_levels,
                "gtop_reactions": gtop_reactions,
            }.items()
            if value is not None
        }
        if (
            not mapping_inputs
            and pathways is None
            and relations is None
            and complex_pathways is None
            and ewas_pathways is None
            and pathway_gene_sets is None
        ):
            raise ValueError("At least one Reactome input file must be provided")
        normalized_release_version = _normalize_release_version(release_version)
        normalized_mapping_inputs: dict[str, Path] = {}
        for argument_name, file_path in mapping_inputs.items():
            spec = MAPPING_ROLE_BY_ARGUMENT[argument_name]
            normalized_mapping_inputs[argument_name] = _validate_reactome_file(
                file_path,
                label=f"Reactome {spec.role} file",
                forbidden_names=MAPPING_OFFICIAL_FILENAMES - {spec.filename},
            )
        file_pathways = (
            None
            if pathways is None
            else _validate_reactome_file(
                pathways,
                label="Reactome pathways file",
                forbidden_names=MAPPING_OFFICIAL_FILENAMES
                | {"ReactomePathwaysRelation.txt"},
            )
        )
        file_relations = (
            None
            if relations is None
            else _validate_reactome_file(
                relations,
                label="Reactome pathway relations file",
                forbidden_names=MAPPING_OFFICIAL_FILENAMES | {"ReactomePathways.txt"},
            )
        )
        file_complex_pathways = (
            None
            if complex_pathways is None
            else _validate_reactome_file(
                complex_pathways,
                label="Reactome complex pathway file",
            )
        )
        file_ewas_pathways = (
            None
            if ewas_pathways is None
            else _validate_reactome_file(
                ewas_pathways,
                label="Reactome EWAS pathway file",
            )
        )
        file_pathway_gene_sets = (
            None
            if pathway_gene_sets is None
            else _validate_reactome_file(
                pathway_gene_sets,
                label="Reactome pathway gene-set archive",
            )
        )
        result = cls(
            snapshot=_ReactomeSnapshot(
                file_uniprot2reactome=normalized_mapping_inputs.get("uniprot_mapping"),
                file_uniprot_all_levels=normalized_mapping_inputs.get(
                    "uniprot_all_levels"
                ),
                file_uniprot_reactions=normalized_mapping_inputs.get(
                    "uniprot_reactions"
                ),
                file_ncbi_mapping=normalized_mapping_inputs.get("ncbi_mapping"),
                file_ncbi_all_levels=normalized_mapping_inputs.get("ncbi_all_levels"),
                file_ncbi_reactions=normalized_mapping_inputs.get("ncbi_reactions"),
                file_chebi_mapping=normalized_mapping_inputs.get("chebi_mapping"),
                file_chebi_all_levels=normalized_mapping_inputs.get("chebi_all_levels"),
                file_chebi_reactions=normalized_mapping_inputs.get("chebi_reactions"),
                file_gtop_mapping=normalized_mapping_inputs.get("gtop_mapping"),
                file_gtop_all_levels=normalized_mapping_inputs.get("gtop_all_levels"),
                file_gtop_reactions=normalized_mapping_inputs.get("gtop_reactions"),
                file_complex_pathways=file_complex_pathways,
                file_ewas_pathways=file_ewas_pathways,
                file_pathway_gene_sets=file_pathway_gene_sets,
                file_pathways=file_pathways,
                file_relations=file_relations,
            ),
        )
        result._release_version = normalized_release_version
        result._release_version_source = (
            "caller" if normalized_release_version is not None else None
        )
        return result

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> ReactomeDatabase:
        """Open a validated Reactome publication for domain and SQL access.

        Validation reads metadata and catalog schemas only; it does not recount
        the biological relations. The returned handle is pinned to the exact
        file that passed validation.

        Args:
            path: A bioextract Reactome metadata-v2 DuckDB publication.

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
            >>> db.select_ids(["P04637"]).mappings().collect().height > 0  # doctest: +SKIP
            True
        """
        publication_path = Path(path).absolute()
        identity_before = _file_identity(publication_path)
        try:
            tables, release_version, release_version_source = (
                _validate_reactome_publication(publication_path)
            )
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
        result._release_version = release_version
        result._release_version_source = release_version_source
        return result

    @property
    def release_version(self) -> str | None:
        """Return the caller-declared Reactome snapshot identity.

        Examples:
            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
            ...     release_version="96",
            ... )
            >>> db.release_version
            '96'
        """
        return self._release_version

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
            ...     .frames["uniprot_pathway_lowest_level"]
            ...     .collect()["species"]
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
        result._release_version = self._release_version
        result._release_version_source = self._release_version_source
        result._publication_path = self._publication_path
        result._publication_identity = self._publication_identity
        result._publication_tables = self._publication_tables
        return result

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: str = "uniprot",
        target: str = "pathway",
        pathway_level: str | None = None,
    ) -> ReactomeSelection:
        """Create a single-query selection for one Reactome capability.

        Args:
            ids: Caller-supplied identifiers in the selected namespace.
                UniProt pipe values and namespace-specific numeric forms are
                normalized by the shared input normalizer.
            namespace: ``"uniprot"``, ``"ncbi"``, ``"chebi"``, or ``"gtop"``.
            target: ``"pathway"`` or ``"reaction"``.
            pathway_level: ``"lowest_level"`` or ``"all_levels"`` for pathway
                targets; reaction targets require ``None``.

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
            >>> selection.mappings().select("input_id").collect().unique().to_series().to_list()
            ['P04637']
            >>> selection.unmatched_ids().collect().to_dicts()
            [{'input_id': 'MISSING'}]
        """
        namespace, target, pathway_level = _validate_selection_dimensions(
            namespace=namespace,
            target=target,
            pathway_level=pathway_level,
        )
        self._assert_publication_current()
        df_input_ids = _create_namespace_input_frame(ids, namespace)
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
            namespace=namespace,
            target=target,
            pathway_level=pathway_level,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: str = "uniprot",
        target: str = "pathway",
        pathway_level: str | None = None,
    ) -> ReactomeSelection:
        """Create a grouped selection for one Reactome capability.

        Args:
            ids_by_group: Mapping from group label to caller-supplied identifiers.
            namespace: ``"uniprot"``, ``"ncbi"``, ``"chebi"``, or ``"gtop"``.
            target: ``"pathway"`` or ``"reaction"``.
            pathway_level: ``"lowest_level"`` or ``"all_levels"`` for pathway
                targets; reaction targets require ``None``.

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
            ...     selection.mappings()
            ...     .select("group_id", "input_id")
            ...     .unique()
            ...     .collect().to_dicts()
            ... )
            [{'group_id': 'tumor', 'input_id': 'P04637'}]
            >>> selection.unmatched_ids().collect().to_dicts()
            [{'group_id': 'control', 'input_id': 'MISSING'}]
        """
        namespace, target, pathway_level = _validate_selection_dimensions(
            namespace=namespace,
            target=target,
            pathway_level=pathway_level,
        )
        self._assert_publication_current()
        grp_in_frames = _create_namespace_group_input_frames(ids_by_group, namespace)
        return ReactomeSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            _df_group_membership=grp_in_frames.df_group_membership,
            namespace=namespace,
            target=target,
            pathway_level=pathway_level,
        )

    def pathway_mappings(
        self,
        *,
        namespace: str = "uniprot",
        pathway_level: str = "lowest_level",
    ) -> pl.LazyFrame:
        """Return canonical pathway mapping rows for one explicit level.

        The default is the official lowest-level UniProt relation. The
        hierarchy-expanded relation is never substituted or unioned implicitly.

        Examples:
            Request hierarchy-expanded UniProt mappings explicitly:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_all_levels="data/reactome/UniProt2Reactome_All_Levels.txt"
            ... )
            >>> db.pathway_mappings(pathway_level="all_levels")  # doctest: +SKIP
            <LazyFrame ...>
        """
        namespace, _, resolved_pathway_level = _validate_selection_dimensions(
            namespace=namespace,
            target="pathway",
            pathway_level=pathway_level,
        )
        assert resolved_pathway_level is not None
        snapshot = copy.copy(self)
        spec = MAPPING_ROLE_BY_DIMENSIONS[
            (namespace, "pathway", resolved_pathway_level)
        ]
        schema = dict.fromkeys(spec.public_columns, pl.String)
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation="pathway_mappings",
                    namespace=namespace,
                    target="pathway",
                    pathway_level=resolved_pathway_level,
                    schema=schema,
                    request=request,
                ),
            )
        snapshot._require_mapping_capability_for(
            namespace, "pathway", resolved_pathway_level
        )
        file_mapping = snapshot._mapping_source_path_for(
            namespace, "pathway", resolved_pathway_level
        )
        assert file_mapping is not None
        lf_mapping = scan_mapping_role_frame(file_mapping, spec)
        if snapshot.species is not None:
            lf_mapping = lf_mapping.filter(pl.col("species") == snapshot.species)
        return lf_mapping.select(list(schema)).unique().sort(list(schema))

    def reaction_mappings(self, *, namespace: str = "uniprot") -> pl.LazyFrame:
        """Return identifier-to-Reactome-reaction evidence rows lazily.

        Reaction mappings describe membership in a Reactome reaction event. They
        do not imply reaction participants, direction, catalysts, or topology.

        Examples:
            >>> db.reaction_mappings().collect()  # doctest: +SKIP
            shape: (..., 6)
        """
        namespace, _, pathway_level = _validate_selection_dimensions(
            namespace=namespace,
            target="reaction",
            pathway_level=None,
        )
        assert pathway_level is None
        snapshot = copy.copy(self)
        spec = MAPPING_ROLE_BY_DIMENSIONS[(namespace, "reaction", None)]
        schema = dict.fromkeys(spec.public_columns, pl.String)
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation="reaction_mappings",
                    namespace=namespace,
                    target="reaction",
                    pathway_level=None,
                    schema=schema,
                    request=request,
                ),
            )
        snapshot._require_mapping_capability_for(namespace, "reaction", None)
        file_mapping = snapshot._mapping_source_path_for(namespace, "reaction", None)
        assert file_mapping is not None
        lf_mapping = scan_mapping_role_frame(file_mapping, spec)
        if snapshot.species is not None:
            lf_mapping = lf_mapping.filter(pl.col("species") == snapshot.species)
        return lf_mapping.select(list(schema)).unique().sort(list(schema))

    def pathway_genes(self, *, pathway_level: str = "lowest_level") -> pl.LazyFrame:
        """Return distinct Reactome pathway-to-UniProt edges lazily.

        ``pathway_level`` is explicit because all-level enrichment has a
        different biological universe from the preserved lowest-level default.

        Examples:
            >>> db.pathway_genes().collect()  # doctest: +SKIP
            shape: (..., 2)
            >>> db.pathway_genes(pathway_level="all_levels")  # doctest: +SKIP
            shape: (..., 2)
        """
        mapping = self.pathway_mappings(pathway_level=pathway_level)
        schema = dict.fromkeys(["reactome_pathway_id", "uniprot_id"], pl.String)
        return mapping.select(list(schema)).unique().sort(list(schema))

    def _eager_pathway_genes(self, pathway_level: str = "lowest_level") -> pl.DataFrame:
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
            >>> db.pathway_genes().collect().head(2).to_dicts()
            [{'reactome_pathway_id': 'R-HSA-6798695', 'uniprot_id': 'P04637'}, {'reactome_pathway_id': 'R-HSA-6798695', 'uniprot_id': 'Q9Y243'}]
        """
        self._assert_publication_current()
        return extract_term2gene_frame(self._mapping_frame(pathway_level))

    def pathway_names(self) -> pl.LazyFrame:
        """Return Reactome pathway display names and species lazily.

        Examples:
            >>> db.pathway_names().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        snapshot = copy.copy(self)
        schema = dict.fromkeys(
            ["reactome_pathway_id", "pathway_name", "species"],
            pl.String,
        )
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation="pathway_names",
                    schema=schema,
                    request=request,
                ),
            )
        if not snapshot._has_pathway():
            snapshot._raise_missing_capability(
                "Cannot extract Reactome pathways without pathways file",
                "Reactome publication does not contain pathway metadata",
            )
        assert snapshot.snapshot.file_pathways is not None
        lf_pathway = scan_pathway_frame(snapshot.snapshot.file_pathways)
        if snapshot.species is not None:
            lf_pathway = lf_pathway.filter(pl.col("species") == snapshot.species)
        return lf_pathway.select(list(schema)).unique().sort(list(schema))

    def _eager_pathway_names(self) -> pl.DataFrame:
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
            ...     db.pathway_names()
            ...     .select("reactome_pathway_id", "pathway_name")
            ...     .collect()
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'reactome_pathway_id': 'R-HSA-1640170', 'pathway_name': 'Cell Cycle'}]
        """
        self._assert_publication_current()
        if self._df_term2name is None:
            self._df_term2name = extract_term2name_frame(self._pathway_frame())
        return self._df_term2name

    def pathway_relations(self) -> pl.LazyFrame:
        """Return Reactome parent-child pathway edges lazily.

        Examples:
            >>> db.pathway_relations().collect()  # doctest: +SKIP
            shape: (..., 2)
        """
        snapshot = copy.copy(self)
        schema = dict.fromkeys(
            ["parent_reactome_pathway_id", "child_reactome_pathway_id"],
            pl.String,
        )
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation="pathway_relations",
                    schema=schema,
                    request=request,
                ),
            )
        if not snapshot._has_relation():
            snapshot._raise_missing_capability(
                "Cannot read Reactome relations without relations file",
                "Reactome publication does not contain pathway relations",
            )
        assert snapshot.snapshot.file_relations is not None
        lf_relations = scan_relation_frame(snapshot.snapshot.file_relations)
        if snapshot.species is not None:
            if not snapshot._has_pathway():
                snapshot._raise_missing_capability(
                    "Reactome species-scoped relation filtering requires pathways file",
                    "Reactome species-scoped relation filtering requires pathway metadata",
                )
            assert snapshot.snapshot.file_pathways is not None
            lf_pathway_ids = (
                scan_pathway_frame(snapshot.snapshot.file_pathways)
                .filter(pl.col("species") == snapshot.species)
                .select("reactome_pathway_id")
                .unique()
            )
            lf_relations = lf_relations.join(
                lf_pathway_ids.rename(
                    {"reactome_pathway_id": "parent_reactome_pathway_id"}
                ),
                on="parent_reactome_pathway_id",
                how="inner",
            ).join(
                lf_pathway_ids.rename(
                    {"reactome_pathway_id": "child_reactome_pathway_id"}
                ),
                on="child_reactome_pathway_id",
                how="inner",
            )
        return lf_relations.select(list(schema)).unique().sort(list(schema))

    def _eager_pathway_relations(self) -> pl.DataFrame:
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
            >>> db.pathway_relations().collect().head(1).to_dicts()
            [{'parent_reactome_pathway_id': 'R-HSA-1640170', 'child_reactome_pathway_id': 'R-HSA-6798695'}]
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
                        "parent_reactome_pathway_id",
                        "child_reactome_pathway_id",
                    )
                )
            else:
                self._df_relations = filter_relation_frame(
                    self._relation_raw_frame(),
                    self._pathway_frame(),
                )
        return self._df_relations

    def complex_pathways(self) -> pl.LazyFrame:
        """Return Reactome Complex-to-human-pathway membership rows.

        Examples:
            >>> db.complex_pathways().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        return self._entity_pathways_relation(COMPLEX_PATHWAY_ROLE)

    def ewas_pathways(self) -> pl.LazyFrame:
        """Return Reactome EntityWithAccessionedSequence pathway rows.

        Examples:
            >>> db.ewas_pathways().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        return self._entity_pathways_relation(EWAS_PATHWAY_ROLE)

    def pathway_gene_sets(self) -> pl.LazyFrame:
        """Return human Reactome GMT pathway/gene-symbol memberships.

        Examples:
            >>> db.pathway_gene_sets().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        schema = dict.fromkeys(GMT_SOURCE_SPEC["public_columns"], pl.String)
        snapshot = copy.copy(self)
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation=PATHWAY_GENE_SET_ROLE,
                    schema=schema,
                    request=request,
                ),
            )
        if snapshot._gmt_source_path() is None:
            snapshot._raise_missing_capability(
                "Cannot extract Reactome pathway gene sets without GMT archive",
                "Reactome publication does not contain pathway gene sets",
            )
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_gmt_source_batches(
                snapshot,
                schema=schema,
                request=request,
            ),
        )

    def _entity_pathways_relation(self, role: str) -> pl.LazyFrame:
        spec = ENTITY_ROLE_BY_ROLE[role]
        schema = dict.fromkeys(spec["public_columns"], pl.String)
        snapshot = copy.copy(self)
        if snapshot._publication_path is not None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_relation_batches(
                    snapshot,
                    relation=role,
                    schema=schema,
                    request=request,
                ),
            )
        file_path = snapshot._entity_source_path(role)
        if file_path is None:
            snapshot._raise_missing_capability(
                f"Cannot extract Reactome {role} without its source file",
                f"Reactome publication does not contain {role}",
            )
            raise AssertionError("missing Reactome entity capability")
        lf = scan_entity_pathway_frame(
            file_path,
            source_columns=spec["source_columns"],
            public_columns=spec["public_columns"],
            context=f"Reactome {role} file",
        )
        if snapshot.species is not None and snapshot.species != "Homo sapiens":
            return lf.filter(pl.lit(False)).select(list(schema))
        return lf.select(list(schema)).unique().sort(list(schema))

    def _gmt_source_path(self) -> Path | None:
        return self.snapshot.file_pathway_gene_sets

    def _pathway_gene_set_frame(self) -> pl.DataFrame:
        if self._df_pathway_gene_sets is None:
            if self._publication_path is not None:
                if PATHWAY_GENE_SET_ROLE not in self._publication_tables:
                    self._raise_missing_capability(
                        "Cannot extract Reactome pathway gene sets without GMT archive",
                        "Reactome publication does not contain pathway gene sets",
                    )
                frame = self._read_publication_table(PATHWAY_GENE_SET_ROLE)
            else:
                file_path = self._gmt_source_path()
                if file_path is None:
                    self._raise_missing_capability(
                        "Cannot extract Reactome pathway gene sets without GMT archive",
                        "Reactome publication does not contain pathway gene sets",
                    )
                    raise AssertionError("missing Reactome GMT capability")
                frame = read_gmt_frame(
                    file_path,
                    public_columns=GMT_SOURCE_SPEC["public_columns"],
                    context="Reactome GMT archive",
                )
            if self.species is not None and self.species != "Homo sapiens":
                frame = frame.head(0)
            self._df_pathway_gene_sets = frame
        return self._df_pathway_gene_sets

    def build_tidy(self) -> TidyDataset:
        """Build a lazy Reactome tidy dataset from the available source files.

        Returns:
            A `TidyDataset` containing only canonical role assets derivable from
            the provided raw files. Species scoping is reflected in every
            applicable frame.

        Examples:
            Build all canonical role frames when every v0.5 source role is
            available:

            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt",
            ...     uniprot_all_levels="data/reactome/UniProt2Reactome_All_Levels.txt",
            ...     pathways="data/reactome/ReactomePathways.txt",
            ...     relations="data/reactome/ReactomePathwaysRelation.txt",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['pathway', 'pathway_relation', 'uniprot_pathway_all_level', 'uniprot_pathway_lowest_level']
        """
        self._assert_publication_current()
        frames: dict[str, pl.DataFrame] = {}
        assets: list[TidyAsset] = []
        for path, kind, frame_name in ASSET_SPECS:
            if not self._has_tidy_role(frame_name):
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
            source_schema_profile=SOURCE_SCHEMA_PROFILE,
            build_id_prefix="reactome-mapping",
            assets=tuple(assets),
            resource_name="reactome",
            release_version=self._release_version,
            release_version_source=self._release_version_source,
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish available Reactome relations as one DuckDB.

        Enrichment projections are omitted because they are derived from the
        canonical role tables.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = ReactomeDatabase.from_files(
            ...     uniprot_mapping="data/reactome/UniProt2Reactome.txt"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(Path(dir_out) / "reactome.duckdb")
            ...     result.tables
            ('uniprot_pathway_lowest_level',)
        """
        if self._publication_path is not None:
            raise CapabilityError(
                "write_duckdb() requires a Reactome source-file handle"
            )
        dataset = self.build_tidy()
        relations = tuple(
            RelationSpec(
                table_name=asset.frame_name,
                frame=dataset.frames[asset.frame_name],
                role=asset.kind,
            )
            for asset in dataset.assets
        )
        sources = tuple(
            SourceFileRecord(
                logical_name=source.logical_name,
                path=source.path,
                media_type=source.media_type,
                bytes=source.bytes,
                sha256=source.sha256,
            )
            for source in self._tidy_sources()
        )
        return write_duckdb_publication(
            relations,
            path,
            resource_name="reactome",
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile=SOURCE_SCHEMA_PROFILE,
            sources=sources,
            release_version=self._release_version,
            release_version_source=self._release_version_source,
            if_exists=if_exists,
            column_mappings=self._column_mappings(),
            validation_issues=self._validation_issues(),
        )

    def _column_mappings(self) -> tuple[tuple[str, str, str, str], ...]:
        mappings: list[tuple[str, str, str, str]] = []
        for role, spec in ENTITY_ROLE_BY_ROLE.items():
            if self._entity_source_path(role) is None:
                continue
            mappings.extend(
                (
                    role,
                    source_column,
                    output_column,
                    ENTITY_COLUMN_MAPPING_REASON,
                )
                for source_column, output_column in zip(
                    spec["source_columns"], spec["public_columns"], strict=True
                )
            )
        return tuple(mappings)

    def _mapping_frame(self, pathway_level: str = "lowest_level") -> pl.DataFrame:
        """Return and cache one normalized UniProt pathway mapping level."""
        return self._mapping_frame_for("uniprot", "pathway", pathway_level)

    def _mapping_frame_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> pl.DataFrame:
        spec = _mapping_spec(namespace, target, pathway_level)
        self._require_mapping_capability_for(namespace, target, pathway_level)
        if spec.role not in self._df_mapping_by_role:
            self._df_mapping_by_role[spec.role] = filter_species_frame(
                self._mapping_raw_frame_for(namespace, target, pathway_level),
                self.species,
            )
        return self._df_mapping_by_role[spec.role]

    def _entity_source_path(self, role: str) -> Path | None:
        spec = ENTITY_ROLE_BY_ROLE[role]
        argument_name = spec["argument_name"]
        return getattr(self.snapshot, _snapshot_field_for_argument(argument_name))

    def _entity_frame(self, role: str) -> pl.DataFrame:
        if role not in ENTITY_ROLE_BY_ROLE:
            raise ValueError(f"Unsupported Reactome entity role: {role}")
        if role not in self._df_entity_by_role:
            spec = ENTITY_ROLE_BY_ROLE[role]
            file_path = self._entity_source_path(role)
            if file_path is None:
                if role in self._publication_tables:
                    frame = self._read_publication_table(role)
                else:
                    self._raise_missing_capability(
                        f"Cannot extract Reactome {role} without its source file",
                        f"Reactome publication does not contain {role}",
                    )
                    raise AssertionError("missing Reactome entity capability")
            else:
                frame = read_entity_pathway_frame(
                    file_path,
                    source_columns=spec["source_columns"],
                    public_columns=spec["public_columns"],
                    context=f"Reactome {role} file",
                )
            if self.species is not None and self.species != "Homo sapiens":
                frame = frame.head(0)
            self._df_entity_by_role[role] = frame
        return self._df_entity_by_role[role]

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

    def _mapping_raw_frame(self, pathway_level: str = "lowest_level") -> pl.DataFrame:
        return self._mapping_raw_frame_for("uniprot", "pathway", pathway_level)

    def _mapping_raw_frame_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> pl.DataFrame:
        spec = _mapping_spec(namespace, target, pathway_level)
        self._require_mapping_capability_for(namespace, target, pathway_level)
        if spec.role not in self._df_mapping_raw_by_role:
            if self._publication_path is None:
                file_mapping = self._mapping_source_path_for(
                    namespace, target, pathway_level
                )
                assert file_mapping is not None
                frame = read_mapping_role_frame(file_mapping, spec)
            else:
                frame = self._read_publication_table(spec.role)
            self._df_mapping_raw_by_role[spec.role] = frame
        return self._df_mapping_raw_by_role[spec.role]

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
                self._df_pathways_raw = self._read_publication_table("pathway")
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
                    "pathway_relation"
                )
        return self._df_relations_raw

    def _mapping_source_path(self, pathway_level: str) -> Path | None:
        return self._mapping_source_path_for("uniprot", "pathway", pathway_level)

    def _mapping_source_path_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> Path | None:
        spec = _mapping_spec(namespace, target, pathway_level)
        return getattr(self.snapshot, _snapshot_field_for_argument(spec.argument_name))

    def _mapping_table(self, pathway_level: str) -> str:
        return self._mapping_table_for("uniprot", "pathway", pathway_level)

    def _mapping_table_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> str:
        return _mapping_spec(namespace, target, pathway_level).role

    def _has_mapping(self, pathway_level: str = "lowest_level") -> bool:
        return self._has_mapping_for("uniprot", "pathway", pathway_level)

    def _has_mapping_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> bool:
        table_name = self._mapping_table_for(namespace, target, pathway_level)
        return self._mapping_source_path_for(
            namespace, target, pathway_level
        ) is not None or (table_name in self._publication_tables)

    def _require_mapping_capability(self, pathway_level: str) -> None:
        self._require_mapping_capability_for("uniprot", "pathway", pathway_level)

    def _require_mapping_capability_for(
        self,
        namespace: str,
        target: str,
        pathway_level: str | None,
    ) -> None:
        if self._has_mapping_for(namespace, target, pathway_level):
            return
        spec = _mapping_spec(namespace, target, pathway_level)
        if target == "pathway" and pathway_level == "all_levels":
            source_message = (
                f"Cannot extract {namespace} all-level Reactome mapping without "
                f"{spec.filename}"
            )
            publication_message = (
                "Reactome publication does not contain all-level "
                f"{namespace} pathway mappings"
            )
        elif target == "pathway":
            source_label = (
                spec.filename[:-4] if spec.filename.endswith(".txt") else spec.filename
            )
            source_message = (
                f"Cannot extract {namespace} Reactome mapping without "
                f"{source_label} file"
            )
            publication_message = (
                "Reactome publication does not contain lowest-level "
                f"{namespace} pathway mappings"
            )
        else:
            source_message = (
                f"Cannot extract {namespace} Reactome reaction mappings without "
                f"{spec.filename}"
            )
            publication_message = (
                f"Reactome publication does not contain {namespace} reaction mappings"
            )
        self._raise_missing_capability(source_message, publication_message)

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

    def _has_tidy_role(self, frame_name: str) -> bool:
        if frame_name in MAPPING_ROLE_BY_ROLE:
            spec = MAPPING_ROLE_BY_ROLE[frame_name]
            return self._has_mapping_for(
                spec.namespace, spec.target, spec.pathway_level
            )
        if frame_name == PATHWAY_ROLE:
            return self._has_pathway()
        if frame_name == RELATION_ROLE:
            return self._has_relation()
        if frame_name in ENTITY_ROLE_BY_ROLE:
            return self._entity_source_path(frame_name) is not None or (
                frame_name in self._publication_tables
            )
        if frame_name == PATHWAY_GENE_SET_ROLE:
            return self._gmt_source_path() is not None or (
                frame_name in self._publication_tables
            )
        raise ValueError(f"Unsupported Reactome tidy role: {frame_name}")

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

    def _read_publication_table(self, table_name: str) -> pl.DataFrame:
        with self.connect() as connection:
            frame = pl.read_database(  # pyright: ignore[reportUnknownMemberType]
                f'SELECT * FROM "{table_name}"', connection
            )
        return frame

    def _validation_issues(self) -> tuple[ValidationIssue, ...]:
        """Return bounded non-fatal publication issues and enforce v0.5 parity."""
        issues: list[ValidationIssue] = []
        if self._has_pathway():
            pathway_ids = self._pathway_frame().select("reactome_pathway_id").unique()
            for spec in MAPPING_ROLE_SPECS:
                if spec.target != "pathway" or not self._has_mapping_for(
                    spec.namespace, spec.target, spec.pathway_level
                ):
                    continue
                role = spec.role
                missing = (
                    self._mapping_frame_for(
                        spec.namespace, spec.target, spec.pathway_level
                    )
                    .select(spec.event_column)
                    .unique()
                    .join(pathway_ids, on=spec.event_column, how="anti")
                    .sort(spec.event_column)
                )
                issues.extend(
                    ValidationIssue(
                        severity="warning",
                        issue_code="missing_pathway_metadata",
                        source_name=role,
                        relation_name=role,
                        identifier_namespace=spec.event_column,
                        identifier_value=str(pathway_id),
                        referenced_relation="pathway",
                        referenced_identifier=str(pathway_id),
                        message=(
                            f"Reactome pathway {pathway_id} is present in {role} "
                            "but absent from pathway metadata"
                        ),
                    )
                    for pathway_id in missing.get_column(spec.event_column).to_list()
                )

            pathway_metadata = (
                self._pathway_raw_frame()
                .select("reactome_pathway_id", "species")
                .unique(subset=["reactome_pathway_id"])
            )
            metadata_by_id = {
                str(row["reactome_pathway_id"]): str(row["species"])
                for row in pathway_metadata.to_dicts()
            }
            for role, _entity_spec in ENTITY_ROLE_BY_ROLE.items():
                if not self._has_tidy_role(role):
                    continue
                entity_frame = self._entity_frame(role)
                for endpoint_column, endpoint_role in (
                    ("reactome_pathway_id", "pathway"),
                    ("top_level_reactome_pathway_id", "top_level_pathway"),
                ):
                    endpoint_values = (
                        entity_frame.select(endpoint_column)
                        .unique()
                        .get_column(endpoint_column)
                        .to_list()
                    )
                    for endpoint_id in endpoint_values:
                        endpoint_text = str(endpoint_id)
                        species = metadata_by_id.get(endpoint_text)
                        if species == "Homo sapiens":
                            continue
                        issues.append(
                            ValidationIssue(
                                severity="warning",
                                issue_code="missing_pathway_metadata",
                                source_name=role,
                                relation_name=role,
                                identifier_namespace=endpoint_role,
                                identifier_value=endpoint_text,
                                referenced_relation=PATHWAY_ROLE,
                                referenced_identifier=endpoint_text,
                                message=(
                                    f"Reactome {endpoint_role} {endpoint_text} "
                                    f"is absent from human pathway metadata for {role}"
                                ),
                            )
                        )

            if self._has_tidy_role(PATHWAY_GENE_SET_ROLE):
                gmt_pathways = (
                    self._pathway_gene_set_frame()
                    .select("reactome_pathway_id")
                    .unique()
                    .get_column("reactome_pathway_id")
                    .to_list()
                )
                for pathway_id in gmt_pathways:
                    pathway_text = str(pathway_id)
                    if metadata_by_id.get(pathway_text) == "Homo sapiens":
                        continue
                    issues.append(
                        ValidationIssue(
                            severity="warning",
                            issue_code="missing_pathway_metadata",
                            source_name=PATHWAY_GENE_SET_ROLE,
                            relation_name=PATHWAY_GENE_SET_ROLE,
                            identifier_namespace="reactome_pathway_id",
                            identifier_value=pathway_text,
                            referenced_relation=PATHWAY_ROLE,
                            referenced_identifier=pathway_text,
                            message=(
                                f"Reactome pathway {pathway_text} is absent from "
                                "human pathway metadata for GMT gene sets"
                            ),
                        )
                    )

        if self._has_relation():
            for namespace in ("uniprot", "ncbi", "chebi", "gtop"):
                if self._has_mapping_for(
                    namespace, "pathway", "lowest_level"
                ) and self._has_mapping_for(namespace, "pathway", "all_levels"):
                    self._validate_mapping_level_closure(namespace)
            if any(self._has_tidy_role(role) for role in ENTITY_ROLE_BY_ROLE):
                self._validate_entity_pathway_hierarchy()
        return tuple(issues)

    def _validate_uniprot_level_closure(self) -> None:
        self._validate_mapping_level_closure("uniprot")

    def _validate_mapping_level_closure(self, namespace: str) -> None:
        low_spec = _mapping_spec(namespace, "pathway", "lowest_level")
        key_columns = [
            low_spec.source_column,
            "reactome_pathway_id",
            "evidence_code",
            "species",
        ]
        relations = self._relation_raw_frame().unique()
        closure = (
            self._mapping_frame_for(namespace, "pathway", "lowest_level")
            .select(key_columns)
            .unique()
        )
        frontier = closure
        while frontier.height:
            expanded = (
                frontier.join(
                    relations,
                    left_on="reactome_pathway_id",
                    right_on="child_reactome_pathway_id",
                    how="inner",
                )
                .select(
                    low_spec.source_column,
                    pl.col("parent_reactome_pathway_id").alias("reactome_pathway_id"),
                    "evidence_code",
                    "species",
                )
                .unique()
                .join(closure, on=key_columns, how="anti")
            )
            if not expanded.height:
                break
            closure = pl.concat([closure, expanded], how="vertical_relaxed").unique()
            frontier = expanded

        official = (
            self._mapping_frame_for(namespace, "pathway", "all_levels")
            .select(key_columns)
            .unique()
        )
        derived_only = closure.join(official, on=key_columns, how="anti").height
        official_only = official.join(closure, on=key_columns, how="anti").height
        if derived_only or official_only:
            raise IntegrityError(
                f"Reactome {namespace} all-level closure mismatch: "
                f"derived_only={derived_only}, official_only={official_only}"
            )

    def _validate_entity_pathway_hierarchy(self) -> None:
        relations = self._relation_raw_frame().unique().to_dicts()
        parents: dict[str, set[str]] = {}
        nodes: set[str] = set()
        for row in relations:
            parent = str(row["parent_reactome_pathway_id"])
            child = str(row["child_reactome_pathway_id"])
            nodes.update((parent, child))
            parents.setdefault(child, set()).add(parent)
            parents.setdefault(parent, set())

        state: dict[str, int] = {}
        for start in nodes:
            if state.get(start, 0) != 0:
                continue
            stack: list[tuple[str, bool]] = [(start, False)]
            while stack:
                node, exiting = stack.pop()
                if exiting:
                    state[node] = 2
                    continue
                current = state.get(node, 0)
                if current == 1:
                    raise IntegrityError("Reactome pathway hierarchy contains a cycle")
                if current == 2:
                    continue
                state[node] = 1
                stack.append((node, True))
                stack.extend((parent, False) for parent in parents.get(node, ()))

        ancestor_cache: dict[str, set[str]] = {}

        def ancestors(pathway_id: str) -> set[str]:
            if pathway_id in ancestor_cache:
                return ancestor_cache[pathway_id]
            found: set[str] = set()
            pending = list(parents.get(pathway_id, ()))
            while pending:
                parent = pending.pop()
                if parent in found:
                    continue
                found.add(parent)
                pending.extend(parents.get(parent, ()))
            ancestor_cache[pathway_id] = found
            return found

        for role in ENTITY_ROLE_BY_ROLE:
            if not self._has_tidy_role(role):
                continue
            frame = self._entity_frame(role)
            for row in (
                frame.select("reactome_pathway_id", "top_level_reactome_pathway_id")
                .unique()
                .to_dicts()
            ):
                pathway_id = str(row["reactome_pathway_id"])
                top_level_id = str(row["top_level_reactome_pathway_id"])
                if top_level_id not in {pathway_id} | ancestors(pathway_id):
                    raise IntegrityError(
                        f"Reactome {role} top-level pathway is not an ancestor: "
                        f"pathway={pathway_id}, top_level={top_level_id}"
                    )

    def _build_tidy_frame(self, frame_name: str) -> pl.DataFrame:
        if frame_name in MAPPING_ROLE_BY_ROLE:
            spec = MAPPING_ROLE_BY_ROLE[frame_name]
            return self._mapping_frame_for(
                spec.namespace, spec.target, spec.pathway_level
            )
        if frame_name == PATHWAY_ROLE:
            return self._pathway_frame()
        if frame_name == RELATION_ROLE:
            return self._eager_pathway_relations()
        if frame_name in ENTITY_ROLE_BY_ROLE:
            return self._entity_frame(frame_name)
        if frame_name == PATHWAY_GENE_SET_ROLE:
            return self._pathway_gene_set_frame()
        raise ValueError(f"Unsupported Reactome tidy frame: {frame_name}")

    def _tidy_sources(self) -> tuple[TidySource, ...]:
        sources: list[TidySource] = []
        for spec in MAPPING_ROLE_SPECS:
            file_path = self._mapping_source_path_for(
                spec.namespace, spec.target, spec.pathway_level
            )
            if file_path is not None:
                sources.append(
                    TidySource(
                        logical_name=spec.role,
                        path=file_path,
                        media_type=MEDIA_TYPE_TSV,
                    )
                )
        for role, _spec in ENTITY_ROLE_BY_ROLE.items():
            file_path = self._entity_source_path(role)
            if file_path is not None:
                sources.append(
                    TidySource(
                        logical_name=role,
                        path=file_path,
                        media_type=MEDIA_TYPE_TSV,
                    )
                )
        gmt_path = self._gmt_source_path()
        if gmt_path is not None:
            sources.append(
                TidySource(
                    logical_name=PATHWAY_GENE_SET_ROLE,
                    path=gmt_path,
                    media_type=MEDIA_TYPE_ZIP,
                )
            )
        if self.snapshot.file_pathways is not None:
            sources.append(
                TidySource(
                    logical_name=PATHWAY_ROLE,
                    path=self.snapshot.file_pathways,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        if self.snapshot.file_relations is not None:
            sources.append(
                TidySource(
                    logical_name=RELATION_ROLE,
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
        ...     selection.mappings()
        ...     .select("input_id", "reactome_pathway_id")
        ...     .collect()
        ...     .head(1)
        ...     .to_dicts()
        ... )
        [{'input_id': 'P04637', 'reactome_pathway_id': 'R-HSA-6798695'}]
    """

    dataset: ReactomeDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    namespace: str = field(default="uniprot", repr=False)
    target: str = field(default="pathway", repr=False)
    pathway_level: str | None = field(default=None, repr=False)

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

    def mappings(self) -> pl.LazyFrame:
        """Return selected Reactome mapping rows lazily.

        Examples:
            >>> selection.mappings().collect()  # doctest: +SKIP
            shape: (..., ...)
        """
        snapshot = copy.copy(self)
        spec = _mapping_spec(
            snapshot.namespace, snapshot.target, snapshot.pathway_level
        )
        columns = (["group_id"] if self._df_group_membership is not None else []) + [
            "input_id",
            *_selection_public_columns(spec),
        ]
        schema = dict.fromkeys(columns, pl.String)
        if snapshot.dataset._publication_path is not None:  # pyright: ignore[reportPrivateUsage]
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_selection_mapping_batches(
                    snapshot,
                    schema=schema,
                    request=request,
                ),
            )
        snapshot.dataset._require_mapping_capability_for(  # pyright: ignore[reportPrivateUsage]
            snapshot.namespace,
            snapshot.target,
            snapshot.pathway_level,
        )
        file_mapping = snapshot.dataset._mapping_source_path_for(  # pyright: ignore[reportPrivateUsage]
            snapshot.namespace, snapshot.target, snapshot.pathway_level
        )
        assert file_mapping is not None
        lookup_ids = snapshot._df_input_ids.get_column("lookup_id").to_list()  # pyright: ignore[reportPrivateUsage]
        lf_mapping = scan_mapping_role_frame(file_mapping, spec).filter(
            pl.col(spec.source_column).is_in(lookup_ids)
        )
        if snapshot.dataset.species is not None:
            lf_mapping = lf_mapping.filter(
                pl.col("species") == snapshot.dataset.species
            )
        lf_mapping = lf_mapping.with_columns(
            pl.col(spec.source_column).alias("lookup_id")
        )
        lf_input = snapshot._df_input_ids.lazy()  # pyright: ignore[reportPrivateUsage]
        lf_mapping = lf_input.join(lf_mapping, on="lookup_id", how="inner")
        if snapshot._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]
            lf_mapping = (
                snapshot._df_group_membership.lazy()  # pyright: ignore[reportPrivateUsage]
                .join(lf_mapping, on="input_id", how="inner")
                .select(list(schema))
            )
        else:
            lf_mapping = lf_mapping.select(list(schema))
        return lf_mapping.unique().sort(list(schema))

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return selected UniProt IDs without a Reactome mapping lazily.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (
            ["group_id", "input_id"]
            if self._df_group_membership is not None
            else ["input_id"]
        )
        if snapshot.dataset._publication_path is not None:  # pyright: ignore[reportPrivateUsage]
            return register_replayable_source(
                schema=dict.fromkeys(columns, pl.String),
                batches=lambda request: _iter_selection_unmatched_batches(
                    snapshot,
                    request=request,
                ),
            )
        snapshot.dataset._require_mapping_capability_for(  # pyright: ignore[reportPrivateUsage]
            snapshot.namespace,
            snapshot.target,
            snapshot.pathway_level,
        )
        file_mapping = snapshot.dataset._mapping_source_path_for(  # pyright: ignore[reportPrivateUsage]
            snapshot.namespace, snapshot.target, snapshot.pathway_level
        )
        assert file_mapping is not None
        spec = _mapping_spec(
            snapshot.namespace, snapshot.target, snapshot.pathway_level
        )
        lookup_ids = snapshot._df_input_ids.get_column("lookup_id").to_list()  # pyright: ignore[reportPrivateUsage]
        lf_mapping = scan_mapping_role_frame(file_mapping, spec).filter(
            pl.col(spec.source_column).is_in(lookup_ids)
        )
        if snapshot.dataset.species is not None:
            lf_mapping = lf_mapping.filter(
                pl.col("species") == snapshot.dataset.species
            )
        lf_mapping = lf_mapping.select(
            pl.col(spec.source_column).alias("lookup_id")
        ).unique()
        input_rows = (
            snapshot._df_group_membership.lazy()  # pyright: ignore[reportPrivateUsage]
            if snapshot._df_group_membership is not None  # pyright: ignore[reportPrivateUsage]
            else snapshot._df_input_ids.lazy()  # pyright: ignore[reportPrivateUsage]
        )
        return (
            input_rows.join(lf_mapping, on="lookup_id", how="anti")
            .select(columns)
            .sort(columns)
        )


def _validate_selection_dimensions(
    *,
    namespace: str,
    target: str,
    pathway_level: str | None,
) -> tuple[str, str, str | None]:
    if namespace not in {"uniprot", "ncbi", "chebi", "gtop"}:
        raise ValueError(f"Unsupported Reactome namespace: {namespace!r}")
    if target not in {"pathway", "reaction"}:
        raise ValueError(f"Unsupported Reactome target: {target!r}")
    if target == "reaction":
        if pathway_level is not None:
            raise ValueError("Reaction selections do not accept pathway_level")
        return namespace, target, None
    resolved_level = "lowest_level" if pathway_level is None else pathway_level
    if resolved_level not in {"lowest_level", "all_levels"}:
        raise ValueError(f"Unsupported Reactome pathway level: {pathway_level!r}")
    return namespace, target, resolved_level


def _mapping_spec(
    namespace: str,
    target: str,
    pathway_level: str | None,
):
    try:
        return MAPPING_ROLE_BY_DIMENSIONS[(namespace, target, pathway_level)]
    except KeyError as error:
        raise ValueError(
            "Unsupported Reactome mapping dimensions: "
            f"namespace={namespace!r}, target={target!r}, "
            f"pathway_level={pathway_level!r}"
        ) from error


def _selection_public_columns(spec: Any) -> tuple[str, ...]:
    """Return the stable selected-output order, distinct from raw order."""
    return (
        spec.source_column,
        spec.event_column,
        spec.name_column,
        "evidence_code",
        "species",
        "reactome_url",
    )


def _snapshot_field_for_argument(argument_name: str) -> str:
    if argument_name == "uniprot_mapping":
        return "file_uniprot2reactome"
    return f"file_{argument_name}"


_RE_CHEBI_ID = re.compile(r"^(?:chebi:)?([0-9]+)$", re.IGNORECASE)
_RE_DECIMAL_ID = re.compile(r"^[0-9]+$")


def _normalize_namespace_identifier(
    value: object,
    namespace: str,
) -> tuple[str, str] | None:
    raw = str(value).strip()
    if not raw:
        return None
    if namespace == "uniprot":
        normalized = normalize_input_id(raw)
        return (normalized, normalized) if normalized else None
    if namespace == "ncbi":
        return raw, raw
    if namespace == "chebi":
        match = _RE_CHEBI_ID.fullmatch(raw)
        if match is None:
            raise ValueError(f"Invalid ChEBI identifier: {raw!r}")
        numeric = str(int(match.group(1)))
        return f"CHEBI:{numeric}", numeric
    if namespace == "gtop":
        if _RE_DECIMAL_ID.fullmatch(raw) is None:
            raise ValueError(f"Invalid GtoP identifier: {raw!r}")
        return raw, raw
    raise ValueError(f"Unsupported Reactome namespace: {namespace!r}")


def _create_namespace_input_frame(
    input_ids: Iterable[str],
    namespace: str,
) -> pl.DataFrame:
    rows = [
        normalized
        for value in input_ids
        if (normalized := _normalize_namespace_identifier(value, namespace)) is not None
    ]
    schema = {"input_id": pl.String, "lookup_id": pl.String}
    if not rows:
        return pl.DataFrame(schema=schema)
    return (
        pl.DataFrame(
            {
                "input_id": [row[0] for row in rows],
                "lookup_id": [row[1] for row in rows],
            },
            schema=schema,
        )
        .unique(subset=["lookup_id"], maintain_order=True)
        .sort("input_id")
    )


def _create_namespace_group_input_frames(
    ids_by_group: Mapping[str, Iterable[str]],
    namespace: str,
) -> GroupInputFrames:
    group_ids = [str(group_id).strip() for group_id in ids_by_group]
    validate_group_ids(group_ids)
    membership_rows: list[tuple[str, str, str]] = []
    for group_id_raw, ids in ids_by_group.items():
        group_id = str(group_id_raw).strip()
        for value in ids:
            normalized = _normalize_namespace_identifier(value, namespace)
            if normalized is not None:
                membership_rows.append((group_id, normalized[0], normalized[1]))
    group_schema = {"group_id": pl.String}
    membership_schema = {
        "group_id": pl.String,
        "input_id": pl.String,
        "lookup_id": pl.String,
    }
    df_groups = pl.DataFrame({"group_id": group_ids}, schema=group_schema).sort(
        "group_id"
    )
    if not membership_rows:
        empty_membership = pl.DataFrame(schema=membership_schema)
        return GroupInputFrames(
            df_groups=df_groups,
            df_group_membership=empty_membership,
            df_input_ids=pl.DataFrame(
                schema={"input_id": pl.String, "lookup_id": pl.String}
            ),
        )
    df_group_membership = (
        pl.DataFrame(
            {
                "group_id": [row[0] for row in membership_rows],
                "input_id": [row[1] for row in membership_rows],
                "lookup_id": [row[2] for row in membership_rows],
            },
            schema=membership_schema,
        )
        .unique()
        .sort("group_id", "input_id", "lookup_id")
    )
    return GroupInputFrames(
        df_groups=df_groups,
        df_group_membership=df_group_membership,
        df_input_ids=df_group_membership.select("input_id", "lookup_id")
        .unique(subset=["lookup_id"])
        .sort("input_id"),
    )


def _normalize_release_version(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        raise ValueError("release_version must be non-empty when provided")
    return normalized


def _validate_reactome_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
    forbidden_names: Collection[str] = (),
) -> Path:
    file_path = Path(file_path)
    if not file_path.is_file():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    if file_path.name in forbidden_names:
        raise ValueError(f"{label} has the wrong declared role: {file_path.name}")
    return file_path


def _iter_gmt_source_batches(
    database: ReactomeDatabase,
    *,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    database._assert_publication_current()  # pyright: ignore[reportPrivateUsage]
    frame = database._pathway_gene_set_frame()  # pyright: ignore[reportPrivateUsage]
    requested = (
        list(schema)
        if request.columns is None
        else [name for name in request.columns if name in schema]
    )
    if requested:
        frame = frame.select(requested)
    for offset in range(0, frame.height, request.effective_batch_size):
        yield frame.slice(offset, request.effective_batch_size).cast(
            {name: schema[name] for name in frame.columns},
            strict=False,
        )


def _iter_publication_relation_batches(
    database: ReactomeDatabase,
    *,
    relation: str,
    namespace: str = "uniprot",
    target: str = "pathway",
    pathway_level: str | None = "lowest_level",
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    database._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired publication boundary
    species = database.species
    if relation in {"pathway_mappings", "reaction_mappings", "pathway_genes"}:
        if relation == "pathway_genes":
            namespace = "uniprot"
            target = "pathway"
            pathway_level = pathway_level or "lowest_level"
        database._require_mapping_capability_for(  # pyright: ignore[reportPrivateUsage]
            namespace, target, pathway_level
        )
        spec = _mapping_spec(namespace, target, pathway_level)
        table_name = database._mapping_table_for(  # pyright: ignore[reportPrivateUsage]
            namespace, target, pathway_level
        )
        if relation == "pathway_genes":
            query = f'SELECT DISTINCT reactome_pathway_id, {spec.source_column} FROM "{table_name}"'
        else:
            query = (
                f'SELECT DISTINCT {", ".join(spec.public_columns)} FROM "{table_name}"'
            )
        params: list[str] = []
        if species is not None:
            query += " WHERE species = ?"
            params.append(species)
        if relation == "pathway_genes":
            query += " ORDER BY reactome_pathway_id, uniprot_id"
        elif relation == "reaction_mappings":
            query += " ORDER BY " + ", ".join(spec.public_columns)
        else:
            query += " ORDER BY " + ", ".join(spec.public_columns)
    elif relation == PATHWAY_GENE_SET_ROLE:
        if not database._has_tidy_role(PATHWAY_GENE_SET_ROLE):  # pyright: ignore[reportPrivateUsage]
            database._raise_missing_capability(  # pyright: ignore[reportPrivateUsage]
                "Cannot extract Reactome pathway gene sets without GMT archive",
                "Reactome publication does not contain pathway gene sets",
            )
        columns = GMT_SOURCE_SPEC["public_columns"]
        query = f'SELECT DISTINCT {", ".join(columns)} FROM "{PATHWAY_GENE_SET_ROLE}"'
        params = []
        if species is not None and species != "Homo sapiens":
            query += " WHERE FALSE"
        query += " ORDER BY " + ", ".join(columns)
    elif relation in ENTITY_ROLE_BY_ROLE:
        entity_spec = ENTITY_ROLE_BY_ROLE[relation]
        if not database._has_tidy_role(relation):  # pyright: ignore[reportPrivateUsage]
            database._raise_missing_capability(  # pyright: ignore[reportPrivateUsage]
                f"Cannot extract Reactome {relation} without its source file",
                f"Reactome publication does not contain {relation}",
            )
        columns = entity_spec["public_columns"]
        query = f'SELECT DISTINCT {", ".join(columns)} FROM "{relation}"'
        params = []
        if species is not None and species != "Homo sapiens":
            query += " WHERE FALSE"
        query += " ORDER BY " + ", ".join(columns)
    elif relation == "pathway_names":
        if not database._has_pathway():  # pyright: ignore[reportPrivateUsage]
            database._raise_missing_capability(  # pyright: ignore[reportPrivateUsage]
                "Cannot extract Reactome pathways without pathways file",
                "Reactome publication does not contain pathway metadata",
            )
        query = (
            "SELECT DISTINCT reactome_pathway_id, pathway_name, species FROM pathway"
        )
        params = []
        if species is not None:
            query += " WHERE species = ?"
            params.append(species)
        query += " ORDER BY reactome_pathway_id"
    elif relation == "pathway_relations":
        if not database._has_relation():  # pyright: ignore[reportPrivateUsage]
            database._raise_missing_capability(  # pyright: ignore[reportPrivateUsage]
                "Cannot read Reactome relations without relations file",
                "Reactome publication does not contain pathway relations",
            )
        if species is not None and not database._has_pathway():  # pyright: ignore[reportPrivateUsage]
            database._raise_missing_capability(  # pyright: ignore[reportPrivateUsage]
                "Reactome species-scoped relation filtering requires pathways file",
                "Reactome species-scoped relation filtering requires pathway metadata",
            )
        query = "SELECT relation.parent_reactome_pathway_id, relation.child_reactome_pathway_id FROM pathway_relation AS relation"
        params = []
        if species is not None:
            query += (
                " JOIN pathway AS parent_pathway"
                " ON parent_pathway.reactome_pathway_id = relation.parent_reactome_pathway_id"
                " JOIN pathway AS child_pathway"
                " ON child_pathway.reactome_pathway_id = relation.child_reactome_pathway_id"
                " WHERE parent_pathway.species = ? AND child_pathway.species = ?"
            )
            params.extend([species, species])
        query += " ORDER BY relation.parent_reactome_pathway_id, relation.child_reactome_pathway_id"
    else:
        raise KeyError(f"Unknown Reactome relation: {relation}")
    requested = (
        None
        if request.columns is None
        else [name for name in request.columns if name in schema]
    )
    if requested:
        projection = ", ".join(f'"{name}"' for name in requested)
        query = f"SELECT {projection} FROM ({query}) AS _bioextract_relation"
    connection = database.connect()
    reader: Any = None
    try:
        result = connection.execute(query, params)
        reader = _reactome_arrow_reader(result, request.effective_batch_size)
        for record_batch in reader:
            frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
            yield frame.cast(  # type: ignore[reportArgumentType]
                {name: schema[name] for name in frame.columns},
                strict=False,
            )
    finally:
        if reader is not None:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
        connection.close()


def _iter_selection_mapping_batches(
    selection: ReactomeSelection,
    *,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    database = selection.dataset
    database._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired publication boundary
    input_rows = selection._df_input_ids.to_dicts()  # pyright: ignore[reportPrivateUsage]
    if not input_rows:
        return
    connection = database.connect()
    reader: Any = None
    try:
        connection.execute(
            "CREATE TEMP TABLE _reactome_input(input_id VARCHAR, lookup_id VARCHAR PRIMARY KEY)"
        )
        connection.executemany(
            "INSERT INTO _reactome_input VALUES (?, ?)",
            [(str(row["input_id"]), str(row["lookup_id"])) for row in input_rows],
        )
        table_name = database._mapping_table_for(  # pyright: ignore[reportPrivateUsage]
            selection.namespace, selection.target, selection.pathway_level
        )
        spec = _mapping_spec(
            selection.namespace, selection.target, selection.pathway_level
        )
        query = f"""
            SELECT DISTINCT
                input.input_id AS input_id,
                {", ".join(f"mapping.{column}" for column in _selection_public_columns(spec))}
            FROM _reactome_input AS input
            JOIN "{table_name}" AS mapping
              ON mapping.{spec.source_column} = input.lookup_id
        """
        params: list[str] = []
        if database.species is not None:
            query += " WHERE mapping.species = ?"
            params.append(database.species)
        query += " ORDER BY input_id, " + ", ".join(_selection_public_columns(spec))
        requested = (
            None
            if request.columns is None
            else [name for name in request.columns if name in schema]
        )
        if requested:
            projection = ", ".join(f'"{name}"' for name in requested)
            query = f"SELECT {projection} FROM ({query}) AS _bioextract_relation"
        result = connection.execute(query, params)
        reader = _reactome_arrow_reader(result, request.effective_batch_size)
        for record_batch in reader:
            frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
            if selection._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]
                frame = (
                    selection._df_group_membership.join(  # pyright: ignore[reportPrivateUsage]
                        frame,
                        on="input_id",
                        how="inner",
                    )
                    .select(list(schema))
                    .unique()
                    .sort(list(schema))
                )
            yield frame.cast(  # type: ignore[reportArgumentType]
                {name: schema[name] for name in frame.columns},
                strict=False,
            )
    finally:
        if reader is not None:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
        connection.close()


def _iter_selection_unmatched_batches(
    selection: ReactomeSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    database = selection.dataset
    database._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired publication boundary
    input_rows = selection._df_input_ids.to_dicts()  # pyright: ignore[reportPrivateUsage]
    if not input_rows:
        return
    connection = database.connect()
    try:
        connection.execute(
            "CREATE TEMP TABLE _reactome_input(input_id VARCHAR, lookup_id VARCHAR PRIMARY KEY)"
        )
        connection.executemany(
            "INSERT INTO _reactome_input VALUES (?, ?)",
            [(str(row["input_id"]), str(row["lookup_id"])) for row in input_rows],
        )
        table_name = database._mapping_table_for(  # pyright: ignore[reportPrivateUsage]
            selection.namespace, selection.target, selection.pathway_level
        )
        spec = _mapping_spec(
            selection.namespace, selection.target, selection.pathway_level
        )
        query = (
            "SELECT DISTINCT input.lookup_id FROM _reactome_input AS input "
            f'JOIN "{table_name}" AS mapping '
            f"ON mapping.{spec.source_column} = input.lookup_id"
        )
        params: list[str] = []
        if database.species is not None:
            query += " WHERE mapping.species = ?"
            params.append(database.species)
        rows = connection.execute(query, params).fetchall()
    finally:
        connection.close()

    mapped_ids = {str(row[0]) for row in rows}
    input_frame = selection._df_input_ids.filter(  # pyright: ignore[reportPrivateUsage]
        ~pl.col("lookup_id").is_in(mapped_ids)
    )
    if selection._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]
        input_frame = (
            selection._df_group_membership.join(  # pyright: ignore[reportPrivateUsage]
                input_frame,
                on="input_id",
                how="inner",
            )
            .select("group_id", "input_id")
            .sort("group_id", "input_id")
        )
    else:
        input_frame = input_frame.select("input_id")
    requested = (
        None
        if request.columns is None
        else [name for name in request.columns if name in input_frame.columns]
    )
    if requested:
        input_frame = input_frame.select(requested)
    for offset in range(0, input_frame.height, request.effective_batch_size):
        yield input_frame.slice(offset, request.effective_batch_size)


def _reactome_arrow_reader(result: Any, batch_size: int) -> Any:
    to_arrow_reader = getattr(result, "to_arrow_reader", None)
    if to_arrow_reader is not None:
        return to_arrow_reader(batch_size)
    return result.fetch_record_batch(rows_per_batch=batch_size)


_REACTOME_TABLE_CONTRACTS: dict[str, tuple[str, str, tuple[tuple[str, str], ...]]] = {
    spec.role: (
        spec.role,
        "canonical",
        tuple((column, "VARCHAR") for column in spec.public_columns),
    )
    for spec in MAPPING_ROLE_SPECS
}
_REACTOME_TABLE_CONTRACTS.update(
    {
        PATHWAY_ROLE: (
            PATHWAY_ROLE,
            "canonical",
            (
                ("reactome_pathway_id", "VARCHAR"),
                ("pathway_name", "VARCHAR"),
                ("species", "VARCHAR"),
            ),
        ),
        RELATION_ROLE: (
            RELATION_ROLE,
            "canonical",
            (
                ("parent_reactome_pathway_id", "VARCHAR"),
                ("child_reactome_pathway_id", "VARCHAR"),
            ),
        ),
    }
)
for _entity_role, _entity_spec in ENTITY_ROLE_BY_ROLE.items():
    _REACTOME_TABLE_CONTRACTS[_entity_role] = (
        _entity_role,
        "canonical",
        tuple((column, "VARCHAR") for column in _entity_spec["public_columns"]),
    )
_REACTOME_TABLE_CONTRACTS[PATHWAY_GENE_SET_ROLE] = (
    PATHWAY_GENE_SET_ROLE,
    "canonical",
    tuple((column, "VARCHAR") for column in GMT_SOURCE_SPEC["public_columns"]),
)


def _validate_reactome_publication(
    path: Path,
) -> tuple[frozenset[str], str | None, str | None]:
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
            if (
                metadata.get("bioextract.metadata_schema_version")
                != METADATA_SCHEMA_VERSION
            ):
                raise ValueError("Unsupported Reactome metadata schema version")
            validate_duckdb_metadata_v2(connection, metadata)
            if metadata.get("bioextract.resource_name") != "reactome":
                raise ValueError("DuckDB file is not a bioextract Reactome publication")
            if (
                metadata.get("bioextract.source_schema_profile")
                != SOURCE_SCHEMA_PROFILE
            ):
                raise ValueError("Unsupported Reactome source schema profile")
            if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported Reactome resource schema version")
            if "bioextract.scope" in metadata:
                raise ValueError("Reactome publication scope is unsupported")
            release_version = metadata.get("bioextract.release_version")
            release_version_source = metadata.get("bioextract.release_version_source")
            if release_version is not None and release_version_source != "caller":
                raise ValueError("Reactome publication release identity is unsupported")

            source_rows = connection.execute(
                "SELECT logical_name, bytes, media_type FROM _bioextract.source_file"
            ).fetchall()
            source_roles = {str(row[0]) for row in source_rows}
            allowed_roles = {
                contract[0] for contract in _REACTOME_TABLE_CONTRACTS.values()
            }
            allowed_media_types = {
                PATHWAY_GENE_SET_ROLE: MEDIA_TYPE_ZIP,
                **{
                    role: MEDIA_TYPE_TSV
                    for role in _REACTOME_TABLE_CONTRACTS
                    if role != PATHWAY_GENE_SET_ROLE
                },
            }
            if (
                not source_roles
                or len(source_roles) != len(source_rows)
                or not source_roles <= allowed_roles
                or any(
                    (row[1] is not None and int(row[1]) < 0)
                    or str(row[2]) != allowed_media_types.get(str(row[0]))
                    for row in source_rows
                )
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
                (
                    role,
                    source_column,
                    output_column,
                    ENTITY_COLUMN_MAPPING_REASON,
                )
                for role, spec in ENTITY_ROLE_BY_ROLE.items()
                if role in source_roles
                for source_column, output_column in zip(
                    spec["source_columns"], spec["public_columns"], strict=True
                )
            }
            if column_mappings != expected_column_mappings:
                raise ValueError("Reactome column provenance inventory is unsupported")
    except duckdb.Error as error:
        raise ValueError(f"Cannot open Reactome DuckDB publication: {path}") from error
    return frozenset(expected_tables), release_version, release_version_source


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
