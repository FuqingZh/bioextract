from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from enum import StrEnum
from pathlib import Path
from typing import Literal, cast, overload

import duckdb
import polars as pl

from bioextract._publication import (
    METADATA_SCHEMA_VERSION,
    DuckDBWriteResult,
    validate_duckdb_metadata_v2,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError

from .brite.constant import (
    ASSET_SPECS as BRITE_ASSET_SPECS,
)
from .brite.constant import (
    MEDIA_TYPE_JSON,
    SCHEMA_BRITE,
)
from .brite.constant import (
    SCHEMA_VERSION as BRITE_SCHEMA_VERSION,
)
from .brite.tidy import build_tidy_frames as build_brite_tidy_frames
from .mapping.constant import KEGGNamespace
from .mapping.publication import (
    open_mapping_publication,
    write_mapping_publication,
)
from .mapping.query import (
    KeggSelection,
)
from .mapping.query import (
    gene_pathways as mapping_gene_pathways,
)
from .mapping.query import (
    gene_pathways_via_ko as mapping_gene_pathways_via_ko,
)
from .mapping.query import (
    ko_pathways as mapping_ko_pathways,
)
from .mapping.query import (
    relation as mapping_relation,
)
from .mapping.query import (
    validate_namespace as validate_mapping_namespace,
)
from .mapping.source import (
    MappingSnapshot,
    create_directory_snapshot,
    create_files_snapshot,
)
from .metabolic.core import (
    KEGGMetabolicNamespace,
    KEGGMetabolicSelection,
    MetabolicPublication,
    MetabolicSnapshot,
    validate_selection_namespace,
)
from .metabolic.core import (
    evaluate_modules as evaluate_metabolic_modules,
)
from .metabolic.core import (
    from_metabolic_files as create_metabolic_snapshot,
)
from .metabolic.core import (
    open_publication as open_metabolic_publication,
)
from .metabolic.core import (
    write_duckdb as write_metabolic_duckdb,
)

__all__ = ["KEGGDatabase"]


class _KeggSnapshotKind(StrEnum):
    BRITE_JSON = "brite_json"
    MAPPING_FILES = "mapping_files"
    METABOLIC_FILES = "metabolic_files"
    METABOLIC_PUBLICATION = "metabolic_publication"
    BRITE_PUBLICATION = "brite_publication"
    MAPPING_PUBLICATION = "mapping_publication"


@dataclass(frozen=True, slots=True)
class _KeggSnapshot:
    kind: _KeggSnapshotKind
    file_brite_json: Path | None = None
    mapping: MappingSnapshot | None = None
    metabolic: MetabolicSnapshot | None = None


@dataclass(slots=True)
class KEGGDatabase:
    """Path-first access to a local KEGG resource snapshot.

    A handle represents a BRITE hierarchy, an organism mapping, or a metabolic
    source/publication. Metabolic handles stream official flat records into a
    relational DuckDB and expose reaction-centered domain selections.

    Examples:
        Build a BRITE pathway snapshot:

        >>> brite = KEGGDatabase.from_brite_json("data/kegg/tcar00001.json")
        >>> sorted(brite.build_tidy().frames)
        ['pathway']

        Select UniProt IDs from an organism mapping snapshot:

        >>> mapping = KEGGDatabase.from_mapping_files(
        ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
        ...     gene_ko="data/kegg/gene_ko.tsv",
        ...     gene_pathway="data/kegg/gene_pathway.tsv",
        ...     organism_code="hsa",
        ... )
        >>> mapping.select_ids(
        ...     ["P12345"], namespace="uniprot"
        ... ).matches().select("kegg_gene_id").collect().to_series().to_list()
        ['hsa:1']
    """

    snapshot: _KeggSnapshot
    _metabolic_publication: MetabolicPublication | None = field(
        default=None, init=False, repr=False
    )
    _publication_path: Path | None = field(default=None, init=False, repr=False)

    @classmethod
    def from_metabolic_files(
        cls,
        source: os.PathLike[str] | str | None = None,
        *,
        compound_list: os.PathLike[str] | str | None = None,
        compound_entries: (
            os.PathLike[str] | str | Sequence[os.PathLike[str] | str] | None
        ) = None,
        reaction_list: os.PathLike[str] | str | None = None,
        reaction_entries: (
            os.PathLike[str] | str | Sequence[os.PathLike[str] | str] | None
        ) = None,
        enzyme_list: os.PathLike[str] | str | None = None,
        enzyme_entries: (
            os.PathLike[str] | str | Sequence[os.PathLike[str] | str] | None
        ) = None,
        module_list: os.PathLike[str] | str | None = None,
        module_entries: (
            os.PathLike[str] | str | Sequence[os.PathLike[str] | str] | None
        ) = None,
        compound_pubchem: os.PathLike[str] | str | None = None,
        compound_reaction: os.PathLike[str] | str | None = None,
        reaction_enzyme: os.PathLike[str] | str | None = None,
        reaction_ko: os.PathLike[str] | str | None = None,
        reaction_module: os.PathLike[str] | str | None = None,
        reaction_pathway: os.PathLike[str] | str | None = None,
        module_pathway: os.PathLike[str] | str | None = None,
        release_version: str | None = None,
    ) -> KEGGDatabase:
        """Create a partial or complete metabolic handle from local inputs.

        Entry collections may be a directory, one batch, or a sequence of
        batches. Without ``source``, missing roles become absent publication
        capabilities. With ``source``, the directory or zip/tar archive must
        resolve to one complete release after explicit overlays. Every
        explicit non-``None`` role replaces its complete discovered role;
        only ``None`` permits discovery.

        ``release_version`` is exclusively caller-declared. Source paths and
        archive members never supply or validate release identity.

        Raises:
            FileNotFoundError: If a supplied path does not exist.
            ValueError: If no input is supplied, a role resolves to no files,
                the source layout is missing or ambiguous, the merged release
                is incomplete, or two inputs reuse one physical file.

        Examples:
            Replace one role while discovering the rest of a release:

            >>> db = KEGGDatabase.from_metabolic_files(  # doctest: +SKIP
            ...     "kegg/metabolic/2026-07",
            ...     reaction_ko="overrides/reaction_ko.tsv",
            ...     release_version="2026-07",
            ... )
            >>> db.snapshot.metabolic.complete_release  # doctest: +SKIP
            True
        """
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.METABOLIC_FILES,
                metabolic=create_metabolic_snapshot(
                    source=source,
                    compound_list=compound_list,
                    compound_entries=compound_entries,
                    reaction_list=reaction_list,
                    reaction_entries=reaction_entries,
                    enzyme_list=enzyme_list,
                    enzyme_entries=enzyme_entries,
                    module_list=module_list,
                    module_entries=module_entries,
                    compound_pubchem=compound_pubchem,
                    compound_reaction=compound_reaction,
                    reaction_enzyme=reaction_enzyme,
                    reaction_ko=reaction_ko,
                    reaction_module=reaction_module,
                    reaction_pathway=reaction_pathway,
                    module_pathway=module_pathway,
                    release_version=release_version,
                ),
            )
        )

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> KEGGDatabase:
        """Open a validated KEGG publication for domain and read-only access.

        Examples:
            >>> db = KEGGDatabase.from_duckdb("kegg.duckdb")  # doctest: +SKIP
            >>> db.snapshot.kind.value  # doctest: +SKIP
            'metabolic_publication'
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM reaction").fetchone()[0]
        """
        publication_path = Path(path)
        profile = _read_publication_profile(publication_path)
        if profile == "kegg-metabolic-flat-files-v1":
            publication = open_metabolic_publication(publication_path)
            result = cls(
                snapshot=_KeggSnapshot(kind=_KeggSnapshotKind.METABOLIC_PUBLICATION)
            )
            result._metabolic_publication = publication
            result._publication_path = publication_path
            return result
        if profile == "kegg-organism-mapping-files-v2":
            result = cls(
                snapshot=_KeggSnapshot(
                    kind=_KeggSnapshotKind.MAPPING_PUBLICATION,
                    mapping=open_mapping_publication(publication_path),
                ),
            )
            result._publication_path = publication_path
            return result
        kind = _validate_tidy_publication(publication_path, profile=profile)
        result = cls(snapshot=_KeggSnapshot(kind=kind))
        result._publication_path = publication_path
        return result

    @classmethod
    def from_brite_json(
        cls,
        path: os.PathLike[str] | str,
    ) -> KEGGDatabase:
        """Create a dataset handle from a local KEGG BRITE JSON file.

        Args:
            path: KEGG BRITE hierarchy in JSON form.

        Returns:
            A BRITE-mode handle that can build or write the pathway tidy asset.

        Raises:
            FileNotFoundError: If the JSON file does not exist.

        Examples:
            Open a compact BRITE hierarchy and read its first pathway entry:

            >>> db = KEGGDatabase.from_brite_json("data/kegg/tcar00001.json")
            >>> db.build_tidy().frames["pathway"].select(
            ...     "pathway_level3_kegg_id", "entry_id"
            ... ).head(1).collect().to_dicts()
            [{'pathway_level3_kegg_id': 'tcar00010', 'entry_id': 'U0034_04525'}]
        """
        file_brite_json = _validate_file(
            path,
            label="KEGG BRITE JSON file",
        )
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.BRITE_JSON,
                file_brite_json=file_brite_json,
            ),
        )

    @classmethod
    def from_mapping_directory(
        cls,
        source: os.PathLike[str] | str,
        *,
        organism_list: os.PathLike[str] | str | None = None,
        ko_pathway: os.PathLike[str] | str | None = None,
        release_version: str | None = None,
    ) -> KEGGDatabase:
        """Create a lazy multi-organism mapping handle from one direct root.

        Construction validates the root and optional global overlays but does
        not enumerate organism directories or open biological files. Each
        lazy execution observes the organism directories then available below
        ``source``.

        Examples:
            >>> db = KEGGDatabase.from_mapping_directory("data/kegg/mapping")
            >>> isinstance(db.organisms(), pl.LazyFrame)
            True
        """
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.MAPPING_FILES,
                mapping=create_directory_snapshot(
                    source,
                    organism_list=organism_list,
                    ko_pathway=ko_pathway,
                    release_version=release_version,
                ),
            )
        )

    @classmethod
    def from_mapping_files(
        cls,
        source: os.PathLike[str] | str | None = None,
        *,
        organism_code: str,
        gene_list: os.PathLike[str] | str | None = None,
        uniprot_conversion: os.PathLike[str] | str | None = None,
        ncbi_gene_conversion: os.PathLike[str] | str | None = None,
        gene_ko: os.PathLike[str] | str | None = None,
        gene_pathway: os.PathLike[str] | str | None = None,
        organism_list: os.PathLike[str] | str | None = None,
        ko_pathway: os.PathLike[str] | str | None = None,
        release_version: str | None = None,
    ) -> KEGGDatabase:
        """Create a partial one-organism mapping handle from a direct root,
        explicit role files, or both.

        Explicit files replace conventional children of ``source``. Missing
        roles remain unavailable capabilities; at least one organism role is
        required. Global roles are never discovered from a parent directory.

        Examples:
            >>> db = KEGGDatabase.from_mapping_files(
            ...     "data/kegg/hsa", organism_code="hsa"
            ... )
            >>> isinstance(db.gene_annotations(), pl.LazyFrame)
            True
        """
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.MAPPING_FILES,
                mapping=create_files_snapshot(
                    source,
                    organism_code=organism_code,
                    gene_list=gene_list,
                    uniprot_conversion=uniprot_conversion,
                    ncbi_gene_conversion=ncbi_gene_conversion,
                    gene_ko=gene_ko,
                    gene_pathway=gene_pathway,
                    organism_list=organism_list,
                    ko_pathway=ko_pathway,
                    release_version=release_version,
                ),
            ),
        )

    def with_organisms(self, organism_codes: Sequence[str]) -> KEGGDatabase:
        """Return a mapping handle physically scoped to selected organisms.

        Examples:
            >>> scoped = db.with_organisms(["hsa"])  # doctest: +SKIP
            >>> scoped.organisms().collect()["organism_code"].to_list()  # doctest: +SKIP
            ['hsa']
        """
        mapping = self._require_mapping_snapshot("scope KEGG organisms")
        return replace(
            self,
            snapshot=replace(
                self.snapshot,
                mapping=mapping.with_organisms(organism_codes),
            ),
        )

    def organisms(self) -> pl.LazyFrame:
        """Return organism members and optional official metadata lazily.

        Examples:
            >>> db.organisms().select("organism_code").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        return mapping_relation(
            self._require_mapping_snapshot("read KEGG organisms"), "organism"
        )

    def gene_annotations(self) -> pl.LazyFrame:
        """Return one nested aggregate row per composite KEGG gene lazily.

        Examples:
            >>> db.gene_annotations().select("kegg_gene_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        return mapping_relation(
            self._require_mapping_snapshot("read KEGG gene annotations"),
            "gene_annotation",
        )

    def ko_annotations(self) -> pl.LazyFrame:
        """Return the global KO universe and optional pathway mappings lazily.

        Examples:
            >>> db.ko_annotations().select("ko_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        return mapping_relation(
            self._require_mapping_snapshot("read KEGG KO annotations"),
            "ko_annotation",
        )

    def gene_pathways(self) -> pl.LazyFrame:
        """Return direct gene-to-pathway observations as nested lists.

        Examples:
            >>> db.gene_pathways().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        return mapping_gene_pathways(
            self._require_mapping_snapshot("read KEGG gene pathways")
        )

    def ko_pathways(self) -> pl.LazyFrame:
        """Return direct KO-to-pathway observations as nested lists.

        Examples:
            >>> db.ko_pathways().collect()  # doctest: +SKIP
            shape: (..., 2)
        """
        return mapping_ko_pathways(
            self._require_mapping_snapshot("read KEGG KO pathways")
        )

    def gene_pathways_via_ko(self) -> pl.LazyFrame:
        """Return auditable gene-to-KO-to-pathway traversal results lazily.

        Examples:
            >>> db.gene_pathways_via_ko().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        return mapping_gene_pathways_via_ko(
            self._require_mapping_snapshot("traverse KEGG pathways via KO")
        )

    @overload
    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KEGGMetabolicSelection: ...

    @overload
    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: KEGGNamespace,
        include_obsolete: bool = False,
    ) -> KeggSelection: ...

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: KEGGNamespace | KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KeggSelection | KEGGMetabolicSelection:
        """Create a KEGG mapping selection for one set of input IDs.

        Args:
            ids: Identifiers in the declared mapping or metabolic namespace.
            namespace: Mapping namespaces are ``uniprot``, ``ncbi_gene``, and
                ``kegg_gene``. Metabolic namespaces are validated against the
                relations actually present in the opened publication.
            include_obsolete: For metabolic EC selection, permit exact
                historical deleted/transferred entries instead of applying
                the default accepted-entry policy.

        Returns:
            A selection that can materialize matched rows and unmapped IDs.

        Raises:
            ValueError: If this is a BRITE snapshot or the namespace is invalid.

        Examples:
            Normalize a pipe-style UniProt ID before matching it:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(
            ...     ["sp|P12345|GENE1_HUMAN"], namespace="uniprot"
            ... )
            >>> selection.matches().select(
            ...     "input_id", "kegg_gene_id"
            ... ).unique().collect().to_dicts()
            [{'input_id': 'P12345', 'kegg_gene_id': 'hsa:1'}]
        """
        if self.snapshot.kind == _KeggSnapshotKind.METABOLIC_PUBLICATION:
            publication = self._require_metabolic_publication()
            metabolic_namespace = cast("KEGGMetabolicNamespace", namespace)
            validate_selection_namespace(publication, metabolic_namespace)
            return KEGGMetabolicSelection.from_ids(
                publication=publication,
                ids=ids,
                namespace=metabolic_namespace,
                include_obsolete=include_obsolete,
            )
        mapping = self._require_mapping_snapshot("select KEGG IDs")
        mapping_namespace = cast("KEGGNamespace", namespace)
        validate_mapping_namespace(mapping_namespace)
        return KeggSelection._from_ids(  # pyright: ignore[reportPrivateUsage]
            mapping,
            ids,
            namespace=mapping_namespace,
        )

    @overload
    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KEGGMetabolicSelection: ...

    @overload
    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: KEGGNamespace,
        include_obsolete: bool = False,
    ) -> KeggSelection: ...

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: KEGGNamespace | KEGGMetabolicNamespace,
        include_obsolete: bool = False,
    ) -> KeggSelection | KEGGMetabolicSelection:
        """Create a KEGG mapping selection for named input-ID groups.

        Args:
            ids_by_group: Mapping from group name to IDs in one shared namespace.
                Group names and IDs are normalized before limits are checked.
            namespace: Shared mapping or metabolic namespace. Metabolic
                namespaces are validated against the opened publication's
                actual relation inventory.
            include_obsolete: Apply the metabolic EC historical-entry policy
                independently within every group.

        Returns:
            A selection whose matched and unmapped outputs retain ``group_id``.

        Raises:
            ValueError: If this is a BRITE snapshot, the namespace or a group
                name is invalid.

        Examples:
            Retain the group name on matched mapping rows:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"]}, namespace="uniprot"
            ... )
            >>> selection.matches().select(
            ...     "group_id", "input_id"
            ... ).unique().collect().to_dicts()
            [{'group_id': 'up', 'input_id': 'P12345'}]
        """
        if self.snapshot.kind == _KeggSnapshotKind.METABOLIC_PUBLICATION:
            publication = self._require_metabolic_publication()
            metabolic_namespace = cast("KEGGMetabolicNamespace", namespace)
            validate_selection_namespace(publication, metabolic_namespace)
            return KEGGMetabolicSelection.from_groups(
                publication=publication,
                ids_by_group=ids_by_group,
                namespace=metabolic_namespace,
                include_obsolete=include_obsolete,
            )
        mapping = self._require_mapping_snapshot("select grouped KEGG IDs")
        mapping_namespace = cast("KEGGNamespace", namespace)
        validate_mapping_namespace(mapping_namespace)
        return KeggSelection._from_groups(  # pyright: ignore[reportPrivateUsage]
            mapping,
            ids_by_group,
            namespace=mapping_namespace,
        )

    def build_tidy(self) -> TidyDataset:
        """Build the lazy tidy dataset defined by the snapshot mode.

        Returns:
            A BRITE dataset containing ``pathway`` or a mapping dataset
            containing ``mapping``. Source paths and the mode-specific schema
            version are retained for embedded publication provenance.

        Examples:
            Build a BRITE dataset:

            >>> brite = KEGGDatabase.from_brite_json("data/kegg/tcar00001.json")
            >>> sorted(brite.build_tidy().frames)
            ['pathway']

        """
        if self.snapshot.kind == _KeggSnapshotKind.BRITE_JSON:
            file_brite_json = self._required_path(self.snapshot.file_brite_json)
            frames = {
                frame_name: frame.lazy()
                for frame_name, frame in build_brite_tidy_frames(
                    file_brite_json
                ).items()
            }
            return TidyDataset(
                frames=frames,
                source=TidySource(
                    logical_name="brite_json",
                    path=file_brite_json,
                    media_type=MEDIA_TYPE_JSON,
                ),
                resource_schema_version=BRITE_SCHEMA_VERSION,
                source_schema_profile="kegg-brite-json-v1",
                build_id_prefix=f"kegg-brite-{file_brite_json.stem}",
                assets=tuple(
                    TidyAsset(path=path, kind=kind, frame_name=frame_name)
                    for path, kind, frame_name in BRITE_ASSET_SPECS
                ),
                resource_name="kegg",
                scope="brite",
            )

        if self.snapshot.kind in {
            _KeggSnapshotKind.MAPPING_FILES,
            _KeggSnapshotKind.MAPPING_PUBLICATION,
        }:
            raise CapabilityError(
                "KEGG mapping build_tidy() was removed; use the native lazy "
                "relations or write_duckdb()"
            )
        raise CapabilityError("build_tidy() requires a KEGG BRITE source handle")

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: Literal["fail", "replace"] = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish a KEGG source profile as DuckDB.

        Examples:
            >>> result = db.write_duckdb("kegg.duckdb")  # doctest: +SKIP
            >>> result.path.name  # doctest: +SKIP
            'kegg.duckdb'
        """
        if self.snapshot.kind == _KeggSnapshotKind.BRITE_JSON:
            return self.build_tidy().write_duckdb(
                path,
                if_exists=if_exists,
            )
        if self.snapshot.kind == _KeggSnapshotKind.MAPPING_FILES:
            return write_mapping_publication(
                self._require_mapping_snapshot("publish KEGG mappings"),
                Path(path),
                if_exists=if_exists,
            )
        if self.snapshot.kind != _KeggSnapshotKind.METABOLIC_FILES:
            raise CapabilityError("write_duckdb() requires a KEGG source handle")
        snapshot = self.snapshot.metabolic
        if snapshot is None:
            raise CapabilityError("KEGG metabolic sources are missing")
        return write_metabolic_duckdb(
            snapshot,
            Path(path),
            if_exists=if_exists,
        )

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a new caller-owned native read-only DuckDB connection.

        Examples:
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM reaction").fetchone()[0]
            >>> count >= 0  # doctest: +SKIP
            True
        """
        if self._publication_path is not None:
            return duckdb.connect(str(self._publication_path), read_only=True)
        publication = self._require_metabolic_publication()
        return duckdb.connect(str(publication.path), read_only=True)

    def evaluate_modules(self, ko_ids: Iterable[str]) -> pl.DataFrame:
        """Evaluate exact KEGG module top-level blocks for the supplied KOs.

        The result reports required and satisfied block counts, exact
        completeness, and one-based missing block indexes.

        Examples:
            >>> result = db.evaluate_modules(["K00844", "K12407"])  # doctest: +SKIP
            >>> result.columns  # doctest: +SKIP
            ['module_id', 'required_block_count', 'satisfied_block_count', 'is_complete', 'missing_block_indexes']
        """
        return evaluate_metabolic_modules(self._require_metabolic_publication(), ko_ids)

    def _require_metabolic_publication(self) -> MetabolicPublication:
        if self._metabolic_publication is None:
            raise CapabilityError(
                "KEGG metabolic selection requires a publication-backed handle; "
                "write a DuckDB and reopen it with KEGGDatabase.from_duckdb()"
            )
        return self._metabolic_publication

    def _require_mapping_snapshot(self, action: str) -> MappingSnapshot:
        if self.snapshot.kind not in {
            _KeggSnapshotKind.MAPPING_FILES,
            _KeggSnapshotKind.MAPPING_PUBLICATION,
        }:
            raise ValueError(
                f"Cannot {action} from a KEGG BRITE JSON snapshot or publication"
            )
        if self.snapshot.mapping is None:
            raise CapabilityError("KEGG mapping snapshot is missing")
        return self.snapshot.mapping

    @staticmethod
    def _required_path(path: Path | None) -> Path:
        if path is None:
            raise ValueError("Required KEGG resource path is missing")
        return path


def _validate_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    return file_path


def _read_publication_profile(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with duckdb.connect(str(path), read_only=True) as connection:
            row = connection.execute(
                "SELECT value FROM _bioextract.metadata "
                "WHERE key='bioextract.source_schema_profile'"
            ).fetchone()
    except duckdb.Error as error:
        raise ValueError(f"Cannot open KEGG DuckDB publication: {path}") from error
    if row is None:
        raise ValueError("KEGG publication is missing source schema profile")
    return str(row[0])


def _validate_tidy_publication(
    path: Path,
    *,
    profile: str,
) -> _KeggSnapshotKind:
    profiles = {
        "kegg-brite-json-v1": (
            _KeggSnapshotKind.BRITE_PUBLICATION,
            "brite",
            BRITE_SCHEMA_VERSION,
            {"pathway": "canonical"},
        ),
    }
    expected = profiles.get(profile)
    if expected is None:
        raise ValueError(f"Unsupported KEGG source schema profile: {profile}")
    kind, scope, schema_version, expected_tables = expected
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        if (
            metadata.get("bioextract.metadata_schema_version")
            != METADATA_SCHEMA_VERSION
        ):
            raise ValueError("Unsupported KEGG metadata schema version")
        validate_duckdb_metadata_v2(connection, metadata)
        if (
            metadata.get("bioextract.resource_name") != "kegg"
            or metadata.get("bioextract.scope") != scope
        ):
            raise ValueError(
                f"DuckDB file is not a bioextract KEGG {scope} publication"
            )
        if metadata.get("bioextract.resource_schema_version") != schema_version:
            raise ValueError(f"Unsupported KEGG {scope} resource schema version")
        source_roles = {
            str(row[0])
            for row in connection.execute(
                "SELECT logical_name FROM _bioextract.source_file"
            ).fetchall()
        }
        required_source_roles = {"brite_json"}
        allowed_source_roles = required_source_roles
        if not required_source_roles <= source_roles or not source_roles <= (
            allowed_source_roles
        ):
            raise ValueError(f"KEGG {scope} source role inventory is unsupported")
        recorded_rows = connection.execute(
            "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
        ).fetchall()
        recorded = {str(row[0]): str(row[1]) for row in recorded_rows}
        physical = {
            str(row[0])
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema='main' AND table_type='BASE TABLE'"
            ).fetchall()
        }
        if recorded != expected_tables or physical != set(expected_tables):
            raise ValueError(f"KEGG {scope} table inventory does not match metadata")
        for table_name, _, row_count in recorded_rows:
            expected_columns = list(SCHEMA_BRITE)
            actual_columns = [
                (str(row[1]), str(row[2]))
                for row in connection.execute(
                    f'PRAGMA table_info("{table_name}")'
                ).fetchall()
            ]
            if actual_columns != [(column, "VARCHAR") for column in expected_columns]:
                raise ValueError(f"KEGG {scope} table schema is unsupported")
            actual = connection.execute(
                f'SELECT count(*) FROM "{table_name}"'
            ).fetchone()
            if actual is None or int(actual[0]) != int(row_count):
                raise ValueError(f"KEGG {scope} row-count drift: {table_name}")
        observed_column_mapping = {
            tuple(str(value) for value in row)
            for row in connection.execute(
                "SELECT table_name, source_column, output_column, reason "
                "FROM _bioextract.column_mapping"
            ).fetchall()
        }
        if observed_column_mapping:
            raise ValueError("KEGG BRITE column provenance inventory is unsupported")
    return kind
