from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal, cast, overload

import duckdb
import polars as pl

from bioextract._publication import DuckDBWriteResult, ParquetWriteResult
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError

from .brite.constant import (
    ASSET_SPECS as BRITE_ASSET_SPECS,
)
from .brite.constant import (
    MEDIA_TYPE_JSON,
)
from .brite.constant import (
    SCHEMA_VERSION as BRITE_SCHEMA_VERSION,
)
from .brite.tidy import build_tidy_frames as build_brite_tidy_frames
from .mapping.constant import (
    ASSET_SPECS as MAPPING_ASSET_SPECS,
)
from .mapping.constant import (
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    KEGGNamespace,
)
from .mapping.constant import (
    SCHEMA_VERSION as MAPPING_SCHEMA_VERSION,
)
from .mapping.util import (
    build_mapping_frame,
    extract_mapping_frame,
    extract_unmatched_ids_frame,
    read_conv_ncbi_geneid_frame,
    read_conv_uniprot_frame,
    read_gene_ko_frame,
    read_gene_list_frame,
    read_gene_pathway_frame,
    validate_namespace,
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
    from_metabolic_release as discover_metabolic_snapshot,
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


@dataclass(frozen=True, slots=True)
class _KeggSnapshot:
    kind: _KeggSnapshotKind
    file_brite_json: Path | None = None
    file_conv_uniprot: Path | None = None
    file_gene_ko: Path | None = None
    file_gene_pathway: Path | None = None
    organism_code: str | None = None
    file_gene_list: Path | None = None
    file_conv_ncbi_geneid: Path | None = None
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
        ... ).extract_mapping()["KeggGeneId"].to_list()
        ['hsa:1', 'hsa:1']
    """

    snapshot: _KeggSnapshot
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _metabolic_publication: MetabolicPublication | None = field(
        default=None, init=False, repr=False
    )

    @classmethod
    def from_metabolic_release(
        cls,
        source: os.PathLike[str] | str,
        *,
        release_version: str | None = None,
    ) -> KEGGDatabase:
        """Discover a complete local KEGG metabolic release.

        ``source`` may be the release directory, its ``raw`` directory, or a
        zip/tar archive containing the layout. ``release_version`` is an
        optional caller-declared official identity. Paths and archive names
        never supply or validate it. No network access is performed.

        Examples:
            >>> db = KEGGDatabase.from_metabolic_release(  # doctest: +SKIP
            ...     "kegg/metabolic/2026-07"
            ... )
            >>> db.snapshot.kind.value  # doctest: +SKIP
            'metabolic_files'
        """
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.METABOLIC_FILES,
                metabolic=discover_metabolic_snapshot(
                    source, release_version=release_version
                ),
            )
        )

    @classmethod
    def from_metabolic_files(
        cls,
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
        """Create a partial or complete metabolic handle from explicit roles.

        Entry collections may be a directory, one batch, or a sequence of
        batches. Missing roles become absent publication capabilities.

        Examples:
            >>> db = KEGGDatabase.from_metabolic_files(  # doctest: +SKIP
            ...     reaction_entries="reaction/",
            ...     reaction_ko="reaction_ko.tsv",
            ...     release_version="2026-07",
            ... )
            >>> db.snapshot.kind.value  # doctest: +SKIP
            'metabolic_files'
        """
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.METABOLIC_FILES,
                metabolic=create_metabolic_snapshot(
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
        """Open a validated KEGG metabolic publication for read-only access.

        Examples:
            >>> db = KEGGDatabase.from_duckdb("kegg.duckdb")  # doctest: +SKIP
            >>> db.snapshot.kind.value  # doctest: +SKIP
            'metabolic_publication'
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM reaction").fetchone()[0]
        """
        publication = open_metabolic_publication(Path(path))
        result = cls(
            snapshot=_KeggSnapshot(kind=_KeggSnapshotKind.METABOLIC_PUBLICATION)
        )
        result._metabolic_publication = publication
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
    def from_mapping_files(
        cls,
        *,
        uniprot_conversion: os.PathLike[str] | str,
        gene_ko: os.PathLike[str] | str,
        gene_pathway: os.PathLike[str] | str,
        organism_code: str,
        gene_list: os.PathLike[str] | str | None = None,
        ncbi_gene_conversion: os.PathLike[str] | str | None = None,
    ) -> KEGGDatabase:
        """Create a dataset handle from one organism's KEGG mapping files.

        The three required files are KEGG ``conv``/``link`` responses for
        UniProt IDs, KO IDs, and pathways. The optional files add NCBI Gene IDs
        and gene display metadata without changing the output schema.

        Args:
            uniprot_conversion: KEGG UniProt-to-gene conversion table.
            gene_ko: KEGG gene-to-KO link table.
            gene_pathway: KEGG gene-to-pathway link table.
            organism_code: KEGG organism code expected as the gene-ID prefix.
            gene_list: Optional KEGG gene list with symbol and description.
            ncbi_gene_conversion: Optional NCBI-Gene-to-KEGG conversion table.

        Returns:
            A mapping-mode handle for extraction, selection, and tidy output.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If ``organism_code`` is empty.

        Examples:
            Open one organism's mapping files and read a normalized gene mapping:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> db.extract_mapping().select(
            ...     "KeggGeneId", "UniProtId", "KoId"
            ... ).row(0, named=True)
            {'KeggGeneId': 'hsa:1', 'UniProtId': 'P12345', 'KoId': 'K00001'}
        """
        organism_code = str(organism_code).strip()
        if not organism_code:
            raise ValueError("KEGG organism_code must be non-empty after normalization")

        file_conv_uniprot = _validate_file(
            uniprot_conversion,
            label="KEGG conv_uniprot file",
        )
        file_gene_ko = _validate_file(
            gene_ko,
            label="KEGG gene_ko file",
        )
        file_gene_pathway = _validate_file(
            gene_pathway,
            label="KEGG gene_pathway file",
        )
        file_gene_list = gene_list
        if file_gene_list is not None:
            file_gene_list = _validate_file(
                file_gene_list,
                label="KEGG gene_list file",
            )
        file_conv_ncbi_geneid = ncbi_gene_conversion
        if file_conv_ncbi_geneid is not None:
            file_conv_ncbi_geneid = _validate_file(
                file_conv_ncbi_geneid,
                label="KEGG conv_ncbi_geneid file",
            )

        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.MAPPING_FILES,
                file_conv_uniprot=file_conv_uniprot,
                file_gene_ko=file_gene_ko,
                file_gene_pathway=file_gene_pathway,
                organism_code=organism_code,
                file_gene_list=file_gene_list,
                file_conv_ncbi_geneid=file_conv_ncbi_geneid,
            ),
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the normalized many-to-many organism mapping.

        Returns:
            One row per distinct joined mapping combination across KEGG gene,
            UniProt, NCBI Gene, KO, and pathway IDs. Columns backed by omitted
            optional files remain nullable.

        Raises:
            ValueError: If called for a BRITE snapshot or if input KEGG gene IDs
                do not match the configured organism code.

        Examples:
            Preserve the two pathway memberships of one KEGG gene:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> db.extract_mapping().filter(
            ...     pl.col("KeggGeneId") == "hsa:1"
            ... ).select("UniProtId", "KeggPathwayId").to_dicts()
            [{'UniProtId': 'P12345', 'KeggPathwayId': 'hsa00010'}, {'UniProtId': 'P12345', 'KeggPathwayId': 'hsa01100'}]
        """
        self._require_mapping_snapshot("extract KEGG mapping")
        if self._df_mapping is None:
            self._df_mapping = build_mapping_frame(
                organism_code=self.snapshot.organism_code or "",
                df_conv_uniprot=read_conv_uniprot_frame(
                    self._required_path(self.snapshot.file_conv_uniprot)
                ),
                df_conv_ncbi_geneid=read_conv_ncbi_geneid_frame(
                    self.snapshot.file_conv_ncbi_geneid
                ),
                df_gene_ko=read_gene_ko_frame(
                    self._required_path(self.snapshot.file_gene_ko)
                ),
                df_gene_pathway=read_gene_pathway_frame(
                    self._required_path(self.snapshot.file_gene_pathway)
                ),
                df_gene_list=read_gene_list_frame(self.snapshot.file_gene_list),
            )
        return self._df_mapping

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
            >>> selection.extract_mapping().select(
            ...     "InputId", "KeggGeneId"
            ... ).unique().to_dicts()
            [{'InputId': 'P12345', 'KeggGeneId': 'hsa:1'}]
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
        self._require_mapping_snapshot("select KEGG IDs")
        mapping_namespace = cast("KEGGNamespace", namespace)
        validate_namespace(mapping_namespace)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return KeggSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
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
            A selection whose matched and unmapped outputs retain ``GroupId``.

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
            >>> selection.extract_mapping().select(
            ...     "GroupId", "InputId"
            ... ).unique().to_dicts()
            [{'GroupId': 'up', 'InputId': 'P12345'}]
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
        self._require_mapping_snapshot("select grouped KEGG IDs")
        mapping_namespace = cast("KEGGNamespace", namespace)
        validate_namespace(mapping_namespace)
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return KeggSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            _df_group_membership=grp_in_frames.df_group_membership,
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

            Build an organism mapping dataset:

            >>> mapping = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> sorted(mapping.build_tidy().frames)
            ['mapping']
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
            )

        self._require_mapping_snapshot("build KEGG mapping tidy dataset")
        return TidyDataset(
            frames={"mapping": self.extract_mapping().lazy()},
            source=self._mapping_tidy_sources(),
            resource_schema_version=MAPPING_SCHEMA_VERSION,
            source_schema_profile="kegg-organism-mapping-files-v1",
            build_id_prefix=f"kegg-mapping-{self.snapshot.organism_code}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in MAPPING_ASSET_SPECS
            ),
            resource_name="kegg",
        )

    def write_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Atomically publish the selected KEGG relation as one Parquet file.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = KEGGDatabase.from_brite_json("data/kegg/tcar00001.json")
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_parquet(Path(dir_out) / "kegg.parquet")
            ...     result.resource_name.startswith("kegg-")
            True
        """
        return self.build_tidy().write_parquet(path, if_exists=if_exists)

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: Literal["fail", "replace"] = "fail",
        include_source_hashes: bool = False,
    ) -> DuckDBWriteResult:
        """Atomically publish a KEGG metabolic snapshot as DuckDB.

        Examples:
            >>> result = db.write_duckdb("kegg.duckdb")  # doctest: +SKIP
            >>> result.path.name  # doctest: +SKIP
            'kegg.duckdb'
        """
        if self.snapshot.kind != _KeggSnapshotKind.METABOLIC_FILES:
            raise CapabilityError(
                "write_duckdb() requires a KEGG metabolic source handle"
            )
        snapshot = self.snapshot.metabolic
        if snapshot is None:
            raise CapabilityError("KEGG metabolic sources are missing")
        return write_metabolic_duckdb(
            snapshot,
            Path(path),
            if_exists=if_exists,
            include_source_hashes=include_source_hashes,
        )

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a new caller-owned native read-only DuckDB connection.

        Examples:
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM reaction").fetchone()[0]
            >>> count >= 0  # doctest: +SKIP
            True
        """
        publication = self._require_metabolic_publication()
        return duckdb.connect(str(publication.path), read_only=True)

    def evaluate_modules(self, ko_ids: Iterable[str]) -> pl.DataFrame:
        """Evaluate exact KEGG module top-level blocks for the supplied KOs.

        The result reports required and satisfied block counts, exact
        completeness, and one-based missing block indexes.

        Examples:
            >>> result = db.evaluate_modules(["K00844", "K12407"])  # doctest: +SKIP
            >>> result.columns  # doctest: +SKIP
            ['ModuleId', 'RequiredBlockCount', 'SatisfiedBlockCount', 'IsComplete', 'MissingBlockIndexes']
        """
        return evaluate_metabolic_modules(self._require_metabolic_publication(), ko_ids)

    def _require_metabolic_publication(self) -> MetabolicPublication:
        if self._metabolic_publication is None:
            raise CapabilityError(
                "KEGG metabolic selection requires a publication-backed handle; "
                "write a DuckDB and reopen it with KEGGDatabase.from_duckdb()"
            )
        return self._metabolic_publication

    def _mapping_tidy_sources(self) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                logical_name="uniprot_conversion",
                path=self._required_path(self.snapshot.file_conv_uniprot),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                logical_name="gene_ko",
                path=self._required_path(self.snapshot.file_gene_ko),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                logical_name="gene_pathway",
                path=self._required_path(self.snapshot.file_gene_pathway),
                media_type=MEDIA_TYPE_TSV,
            ),
        ]
        if self.snapshot.file_gene_list is not None:
            sources.append(
                TidySource(
                    logical_name="gene_list",
                    path=self.snapshot.file_gene_list,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        if self.snapshot.file_conv_ncbi_geneid is not None:
            sources.append(
                TidySource(
                    logical_name="ncbi_gene_conversion",
                    path=self.snapshot.file_conv_ncbi_geneid,
                    media_type=MEDIA_TYPE_TSV,
                )
            )
        return tuple(sources)

    def _require_mapping_snapshot(self, action: str) -> None:
        if self.snapshot.kind != _KeggSnapshotKind.MAPPING_FILES:
            raise ValueError(f"Cannot {action} from a KEGG BRITE JSON snapshot")

    @staticmethod
    def _required_path(path: Path | None) -> Path:
        if path is None:
            raise ValueError("Required KEGG resource path is missing")
        return path


@dataclass(slots=True)
class KeggSelection:
    """Deferred single or grouped query against a KEGG mapping snapshot.

    Selections are created by :meth:`KEGGDatabase.select_ids` or
    :meth:`KEGGDatabase.select_groups`. Matched output retains the normalized
    ``InputId`` and its ``InputNamespace``; grouped selections additionally prepend
    ``GroupId``.

    Examples:
        Materialize matched rows and report IDs that did not map:

        >>> db = KEGGDatabase.from_mapping_files(
        ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
        ...     gene_ko="data/kegg/gene_ko.tsv",
        ...     gene_pathway="data/kegg/gene_pathway.tsv",
        ...     organism_code="hsa",
        ... )
        >>> selection = db.select_ids(
        ...     ["P12345", "MISSING"], namespace="uniprot"
        ... )
        >>> selection.extract_mapping()["KeggGeneId"].unique().to_list()
        ['hsa:1']
        >>> selection.extract_unmatched_ids().to_dicts()
        [{'InputId': 'MISSING'}]
    """

    dataset: KEGGDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    namespace: KEGGNamespace
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"]}, namespace="uniprot"
            ... )
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every KEGG mapping row matched by the selected input IDs.

        Examples:
            Materialize KEGG genes matched by one UniProt accession:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(["P12345"], namespace="uniprot")
            >>> selection.extract_mapping()["KeggGeneId"].to_list()
            ['hsa:1', 'hsa:1']
        """
        if self._df_mapping is None:
            mapping = extract_mapping_frame(
                self.dataset.extract_mapping(),
                self._df_input_ids,
                namespace=self.namespace,
                cols_group_id=(),
            )
            if self._df_group_membership is not None:
                columns = ["GroupId", *mapping.columns]
                mapping = (
                    self._df_group_membership.join(
                        mapping,
                        on="InputId",
                        how="inner",
                    )
                    .select(columns)
                    .unique()
                    .sort(columns)
                )
            self._df_mapping = mapping
        return self._df_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no KEGG mapping row.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain a normalized input accession that did not map:

            >>> db = KEGGDatabase.from_mapping_files(
            ...     uniprot_conversion="data/kegg/conv_uniprot.tsv",
            ...     gene_ko="data/kegg/gene_ko.tsv",
            ...     gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(
            ...     ["P12345", "MISSING"], namespace="uniprot"
            ... )
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            mapping = self.extract_mapping()
            if self._df_group_membership is None:
                self._df_unmapped = extract_unmatched_ids_frame(
                    self._df_input_ids,
                    mapping,
                    cols_group_id=(),
                )
            else:
                mapped_input_ids = mapping.select("InputId").unique()
                self._df_unmapped = (
                    self._df_group_membership.join(
                        mapped_input_ids,
                        on="InputId",
                        how="anti",
                    )
                    .select("GroupId", "InputId")
                    .sort("GroupId", "InputId")
                )
        return self._df_unmapped


def _validate_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    return file_path
