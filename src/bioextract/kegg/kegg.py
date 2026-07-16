from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

from .brite.constant import (
    ASSET_SPECS as BRITE_ASSET_SPECS,
    MEDIA_TYPE_JSON,
    SCHEMA_VERSION as BRITE_SCHEMA_VERSION,
)
from .brite.tidy import build_tidy_frames as build_brite_tidy_frames
from .mapping.constant import (
    ASSET_SPECS as MAPPING_ASSET_SPECS,
    KeggInputIdKind,
    MEDIA_TYPE_TSV,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION as MAPPING_SCHEMA_VERSION,
)
from .mapping.util import (
    build_mapping_frame,
    extract_mapping_frame,
    extract_unmapped_input_ids_frame,
    read_conv_ncbi_geneid_frame,
    read_conv_uniprot_frame,
    read_gene_ko_frame,
    read_gene_list_frame,
    read_gene_pathway_frame,
    validate_kind_input_id,
)

__all__ = [
    "KeggDb",
    "KeggResourceLimits",
    "KeggTidyDataset",
]


class _KeggSnapshotKind(StrEnum):
    BRITE_JSON = "brite_json"
    MAPPING_FILES = "mapping_files"


@dataclass(frozen=True, slots=True)
class KeggResourceLimits:
    """Optional fail-fast limits for KEGG resources and selections.

    File limits are measured in bytes and checked when a snapshot handle is
    created. Count limits apply after input IDs and group names are normalized;
    ``None`` disables the corresponding limit.

    Examples:
        Reject an oversized BRITE snapshot before parsing it:

        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_brite = Path(dir_tmp) / "br08901.json"
        ...     _ = file_brite.write_text("{}", encoding="utf-8")
        ...     limits = KeggResourceLimits(file_brite_json_bytes_max=1)
        ...     try:
        ...         KeggDb.from_brite_json(file_brite, limits=limits)
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_brite_json_bytes_max: int | None = None
    file_conv_uniprot_bytes_max: int | None = None
    file_gene_ko_bytes_max: int | None = None
    file_gene_pathway_bytes_max: int | None = None
    file_gene_list_bytes_max: int | None = None
    file_conv_ncbi_geneid_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


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


KeggTidyDataset = TidyDataset


@dataclass(slots=True)
class KeggDb:
    """Path-first access to a local KEGG resource snapshot.

    A handle represents either a BRITE JSON hierarchy or one organism's KEGG
    mapping files. BRITE handles expose the pathway tidy dataset; mapping
    handles expose the canonical gene mapping, ID selections, and mapping tidy
    dataset. Operations from the other snapshot mode fail explicitly instead
    of interpreting files heuristically.

    Examples:
        Build a BRITE pathway snapshot:

        >>> brite = KeggDb.from_brite_json("data/kegg/tcar00001.json")
        >>> sorted(brite.build_tidy().frames)
        ['pathway']

        Select UniProt IDs from an organism mapping snapshot:

        >>> mapping = KeggDb.from_mapping_files(
        ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
        ...     file_gene_ko="data/kegg/gene_ko.tsv",
        ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
        ...     organism_code="hsa",
        ... )
        >>> mapping.select_ids(
        ...     ["P12345"], kind_input_id="uniprot"
        ... ).extract_mapping()["KeggGeneId"].to_list()
        ['hsa:1', 'hsa:1']
    """

    snapshot: _KeggSnapshot
    limits: KeggResourceLimits
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)

    DEFAULT_RESOURCE_LIMITS = KeggResourceLimits()

    @classmethod
    def from_brite_json(
        cls,
        file_brite_json: os.PathLike[str] | str,
        *,
        limits: KeggResourceLimits | None = None,
    ) -> KeggDb:
        """Create a dataset handle from a local KEGG BRITE JSON file.

        Args:
            file_brite_json: Path to a KEGG BRITE hierarchy in JSON form.
            limits: Optional resource policy. When omitted, no finite size or
                selection limits are imposed.

        Returns:
            A BRITE-mode handle that can build or write the pathway tidy asset.

        Raises:
            FileNotFoundError: If the JSON file does not exist.
            ValueError: If the file exceeds the configured byte limit.

        Examples:
            Open a compact BRITE hierarchy and read its first pathway entry:

            >>> db = KeggDb.from_brite_json("data/kegg/tcar00001.json")
            >>> db.build_tidy().frames["pathway"].select(
            ...     "pathway_level3_kegg_id", "entry_id"
            ... ).head(1).collect().to_dicts()
            [{'pathway_level3_kegg_id': 'tcar00010', 'entry_id': 'U0034_04525'}]
        """
        limits_resolved = KeggResourceLimits() if limits is None else limits
        file_brite_json = _validate_file(
            file_brite_json,
            size_max=limits_resolved.file_brite_json_bytes_max,
            label="KEGG BRITE JSON file",
        )
        return cls(
            snapshot=_KeggSnapshot(
                kind=_KeggSnapshotKind.BRITE_JSON,
                file_brite_json=file_brite_json,
            ),
            limits=limits_resolved,
        )

    @classmethod
    def from_mapping_files(
        cls,
        *,
        file_conv_uniprot: os.PathLike[str] | str,
        file_gene_ko: os.PathLike[str] | str,
        file_gene_pathway: os.PathLike[str] | str,
        organism_code: str,
        file_gene_list: os.PathLike[str] | str | None = None,
        file_conv_ncbi_geneid: os.PathLike[str] | str | None = None,
        limits: KeggResourceLimits | None = None,
    ) -> KeggDb:
        """Create a dataset handle from one organism's KEGG mapping files.

        The three required files are KEGG ``conv``/``link`` responses for
        UniProt IDs, KO IDs, and pathways. The optional files add NCBI Gene IDs
        and gene display metadata without changing the output schema.

        Args:
            file_conv_uniprot: KEGG UniProt-to-gene conversion table.
            file_gene_ko: KEGG gene-to-KO link table.
            file_gene_pathway: KEGG gene-to-pathway link table.
            organism_code: KEGG organism code expected as the gene-ID prefix.
            file_gene_list: Optional KEGG gene list with symbol and description.
            file_conv_ncbi_geneid: Optional NCBI-Gene-to-KEGG conversion table.
            limits: Optional resource policy. When omitted, no finite size or
                selection limits are imposed.

        Returns:
            A mapping-mode handle for extraction, selection, and tidy output.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If ``organism_code`` is empty or a file exceeds its
                configured byte limit.

        Examples:
            Open one organism's mapping files and read a normalized gene mapping:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
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

        limits_resolved = KeggResourceLimits() if limits is None else limits
        file_conv_uniprot = _validate_file(
            file_conv_uniprot,
            size_max=limits_resolved.file_conv_uniprot_bytes_max,
            label="KEGG conv_uniprot file",
        )
        file_gene_ko = _validate_file(
            file_gene_ko,
            size_max=limits_resolved.file_gene_ko_bytes_max,
            label="KEGG gene_ko file",
        )
        file_gene_pathway = _validate_file(
            file_gene_pathway,
            size_max=limits_resolved.file_gene_pathway_bytes_max,
            label="KEGG gene_pathway file",
        )
        if file_gene_list is not None:
            file_gene_list = _validate_file(
                file_gene_list,
                size_max=limits_resolved.file_gene_list_bytes_max,
                label="KEGG gene_list file",
            )
        if file_conv_ncbi_geneid is not None:
            file_conv_ncbi_geneid = _validate_file(
                file_conv_ncbi_geneid,
                size_max=limits_resolved.file_conv_ncbi_geneid_bytes_max,
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
            limits=limits_resolved,
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

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
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

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: KeggInputIdKind,
    ) -> KeggSelection:
        """Create a KEGG mapping selection for one set of input IDs.

        Args:
            ids: UniProt, NCBI Gene, or KEGG gene IDs. Empty values are removed,
                duplicates are folded, and pipe-style UniProt IDs are reduced
                to their accession.
            kind_input_id: Namespace used to join the normalized IDs. Supported
                values are ``uniprot``, ``ncbi_geneid``, and ``kegg_gene``.

        Returns:
            A selection that can materialize matched rows and unmapped IDs.

        Raises:
            ValueError: If this is a BRITE snapshot, the namespace is invalid,
                or the normalized input count exceeds its configured limit.

        Examples:
            Normalize a pipe-style UniProt ID before matching it:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(
            ...     ["sp|P12345|GENE1_HUMAN"], kind_input_id="uniprot"
            ... )
            >>> selection.extract_mapping().select(
            ...     "InputId", "KeggGeneId"
            ... ).unique().to_dicts()
            [{'InputId': 'P12345', 'KeggGeneId': 'hsa:1'}]
        """
        self._require_mapping_snapshot("select KEGG IDs")
        validate_kind_input_id(kind_input_id)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return KeggSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            kind_input_id=kind_input_id,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
        *,
        kind_input_id: KeggInputIdKind,
    ) -> KeggSelection:
        """Create a KEGG mapping selection for named input-ID groups.

        Args:
            group_to_ids: Mapping from group name to IDs in one shared namespace.
                Group names and IDs are normalized before limits are checked.
            kind_input_id: Namespace used to join the normalized IDs. Supported
                values are ``uniprot``, ``ncbi_geneid``, and ``kegg_gene``.

        Returns:
            A selection whose matched and unmapped outputs retain ``GroupId``.

        Raises:
            ValueError: If this is a BRITE snapshot, the namespace or a group
                name is invalid, or a configured group/input limit is exceeded.

        Examples:
            Retain the group name on matched mapping rows:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"]}, kind_input_id="uniprot"
            ... )
            >>> selection.extract_mapping().select(
            ...     "GroupId", "InputId"
            ... ).unique().to_dicts()
            [{'GroupId': 'up', 'InputId': 'P12345'}]
        """
        self._require_mapping_snapshot("select grouped KEGG IDs")
        validate_kind_input_id(kind_input_id)
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
        return KeggSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            kind_input_id=kind_input_id,
        )

    def build_tidy(self) -> KeggTidyDataset:
        """Build the lazy tidy dataset defined by the snapshot mode.

        Returns:
            A BRITE dataset containing ``pathway`` or a mapping dataset
            containing ``mapping``. Source paths and the mode-specific schema
            version are retained for optional manifest generation.

        Examples:
            Build a BRITE dataset:

            >>> brite = KeggDb.from_brite_json("data/kegg/tcar00001.json")
            >>> sorted(brite.build_tidy().frames)
            ['pathway']

            Build an organism mapping dataset:

            >>> mapping = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
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
            return KeggTidyDataset(
                frames=frames,
                source=TidySource(path=file_brite_json, media_type=MEDIA_TYPE_JSON),
                schema_version=BRITE_SCHEMA_VERSION,
                build_id_prefix=f"kegg-brite-{file_brite_json.stem}",
                assets=tuple(
                    TidyAsset(path=path, kind=kind, frame_name=frame_name)
                    for path, kind, frame_name in BRITE_ASSET_SPECS
                ),
            )

        self._require_mapping_snapshot("build KEGG mapping tidy dataset")
        return KeggTidyDataset(
            frames={"mapping": self.extract_mapping().lazy()},
            source=self._mapping_tidy_sources(),
            schema_version=MAPPING_SCHEMA_VERSION,
            build_id_prefix=f"kegg-mapping-{self.snapshot.organism_code}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in MAPPING_ASSET_SPECS
            ),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the mode-specific KEGG tidy assets as flat parquet files.

        Args:
            dir_out: Output directory for the declared tidy assets.
            should_write_manifest: Whether to also write ``manifest.json``.
            should_hash_assets: Whether manifest asset records include SHA-256
                checksums. This has no effect unless a manifest is requested.

        Returns:
            A report containing written asset paths and optional manifest data.

        Examples:
            Write the declared BRITE asset:

            >>> db = KeggDb.from_brite_json("data/kegg/tcar00001.json")
            >>> report = db.write_tidy("build/kegg-brite")
            >>> [asset.path for asset in report.assets]
            ['pathway.parquet']
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def _mapping_tidy_sources(self) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                path=self._required_path(self.snapshot.file_conv_uniprot),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                path=self._required_path(self.snapshot.file_gene_ko),
                media_type=MEDIA_TYPE_TSV,
            ),
            TidySource(
                path=self._required_path(self.snapshot.file_gene_pathway),
                media_type=MEDIA_TYPE_TSV,
            ),
        ]
        if self.snapshot.file_gene_list is not None:
            sources.append(
                TidySource(path=self.snapshot.file_gene_list, media_type=MEDIA_TYPE_TSV)
            )
        if self.snapshot.file_conv_ncbi_geneid is not None:
            sources.append(
                TidySource(
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

    Selections are created by :meth:`KeggDb.select_ids` or
    :meth:`KeggDb.select_groups`. Matched output retains the normalized
    ``InputId`` and its ``KindInputId``; grouped selections additionally prepend
    ``GroupId``.

    Examples:
        Materialize matched rows and report IDs that did not map:

        >>> db = KeggDb.from_mapping_files(
        ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
        ...     file_gene_ko="data/kegg/gene_ko.tsv",
        ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
        ...     organism_code="hsa",
        ... )
        >>> selection = db.select_ids(
        ...     ["P12345", "MISSING"], kind_input_id="uniprot"
        ... )
        >>> selection.extract_mapping()["KeggGeneId"].unique().to_list()
        ['hsa:1']
        >>> selection.extract_unmapped_input_ids().to_dicts()
        [{'InputId': 'MISSING'}]
    """

    dataset: KeggDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: KeggInputIdKind
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"]}, kind_input_id="uniprot"
            ... )
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every KEGG mapping row matched by the selected input IDs.

        Examples:
            Materialize KEGG genes matched by one UniProt accession:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(["P12345"], kind_input_id="uniprot")
            >>> selection.extract_mapping()["KeggGeneId"].to_list()
            ['hsa:1', 'hsa:1']
        """
        if self._df_mapping is None:
            self._df_mapping = extract_mapping_frame(
                self.dataset.extract_mapping(),
                self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no KEGG mapping row.

        Grouped selections report an ID as unmapped independently within each
        group and include ``GroupId`` in the result.

        Examples:
            Retain a normalized input accession that did not map:

            >>> db = KeggDb.from_mapping_files(
            ...     file_conv_uniprot="data/kegg/conv_uniprot.tsv",
            ...     file_gene_ko="data/kegg/gene_ko.tsv",
            ...     file_gene_pathway="data/kegg/gene_pathway.tsv",
            ...     organism_code="hsa",
            ... )
            >>> selection = db.select_ids(
            ...     ["P12345", "MISSING"], kind_input_id="uniprot"
            ... )
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


def _validate_file(
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
