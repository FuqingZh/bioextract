from __future__ import annotations

import json
import os
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
    validate_count_limit,
    validate_file_size,
)
from bioextract._tidy import (
    TidyAsset,
    TidyDataset,
    TidyManifest,
    TidyManifestAsset,
    TidyReportAsset,
    TidySource,
    TidyWriteReport,
    calculate_file_sha256,
)

from ._pfam import build_pfam_tidy_dataset
from .constant import (
    ASSET_SPECS,
    MEDIA_TYPE_TSV_GZIP,
    MEDIA_TYPE_XML_GZIP,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
    InterProInputIdKind,
    InterProTidyConfig,
)
from .util import (
    extract_unmapped_input_ids_frame,
    read_interpro_xml_frames,
    read_mapping_frame,
    scan_mapping_frame,
    select_mapping_frame,
    validate_kind_input_id,
)

__all__ = [
    "InterProDb",
    "InterProResourceLimits",
    "InterProSelection",
    "InterProTidyConfig",
    "InterProTidyDataset",
]


@dataclass(frozen=True, slots=True)
class InterProResourceLimits:
    """Optional guards for InterPro sources and selection cardinality.

    Attributes:
        file_protein2ipr_bytes_max: Maximum on-disk mapping-file size in bytes.
            `None` disables the limit.
        file_interpro_xml_bytes_max: Maximum on-disk XML-file size in bytes.
            `None` disables the limit.
        num_input_ids_max: Maximum number of distinct normalized IDs in one
            selection. `None` disables the limit.
        num_groups_max: Maximum number of groups in one grouped selection.
            `None` disables the limit.

    Examples:
        Reject an oversized mapping snapshot before parsing it:

        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_mapping = Path(dir_tmp) / "protein2ipr.dat"
        ...     _ = file_mapping.write_text("P12345\\tIPR000001\\n")
        ...     limits = InterProResourceLimits(file_protein2ipr_bytes_max=1)
        ...     try:
        ...         InterProDb.from_mapping_files(
        ...             file_protein2ipr=file_mapping, limits=limits
        ...         )
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_protein2ipr_bytes_max: int | None = None
    file_interpro_xml_bytes_max: int | None = None
    num_input_ids_max: int | None = None
    num_groups_max: int | None = None


@dataclass(frozen=True, slots=True)
class _InterProSnapshot:
    file_protein2ipr: Path
    file_interpro_xml: Path | None = None


InterProTidyDataset = TidyDataset


@dataclass(slots=True)
class InterProDb:
    """Access one local InterPro protein-mapping snapshot.

    `protein2ipr.dat.gz` supplies the canonical mapping. Optional InterPro XML
    enriches mapping rows with entry type and member-database metadata; without
    XML those fields remain null. The compact Pfam tidy configuration requires
    same-version XML and validates exact `InterProId + PfamId` relationships.

    Examples:
        Read the first domain annotation from a versioned snapshot:

        >>> db = InterProDb.from_mapping_files(
        ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
        ...     file_interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
        ... )
        >>> db.extract_mapping().select(
        ...     "UniProtId", "InterProId", "MemberDb"
        ... ).head(1).rows()
        [('P12345', 'IPR000001', 'PFAM')]
    """

    snapshot: _InterProSnapshot
    limits: InterProResourceLimits = field(default_factory=InterProResourceLimits)
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _frames_xml: dict[str, pl.DataFrame] | None = field(
        default=None, init=False, repr=False
    )

    DEFAULT_RESOURCE_LIMITS = InterProResourceLimits()

    @classmethod
    def from_mapping_files(
        cls,
        *,
        file_protein2ipr: os.PathLike[str] | str,
        file_interpro_xml: os.PathLike[str] | str | None = None,
        limits: InterProResourceLimits | None = None,
    ) -> InterProDb:
        """Create a handle from local InterPro mapping files.

        Args:
            file_protein2ipr: Path to the required `protein2ipr.dat.gz` source.
            file_interpro_xml: Optional `interpro.xml.gz` source for mapping
                enrichment and required Pfam metadata.
            limits: Optional source-size and selection-cardinality guards.

        Returns:
            A lightweight handle that defers parsing until data is requested.

        Raises:
            FileNotFoundError: If a supplied source file does not exist.
            ValueError: If a supplied file exceeds its configured size limit.

        Examples:
            Use same-version XML to recover entry and member-database metadata:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     file_interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> db.extract_mapping().select(
            ...     "InterProId", "InterProType", "MemberDb"
            ... ).head(1).rows()
            [('IPR000001', 'Domain', 'PFAM')]
        """
        limits_resolved = InterProResourceLimits() if limits is None else limits
        file_protein2ipr = _validate_file(
            file_protein2ipr,
            size_max=limits_resolved.file_protein2ipr_bytes_max,
            label="InterPro protein2ipr file",
        )
        if file_interpro_xml is not None:
            file_interpro_xml = _validate_file(
                file_interpro_xml,
                size_max=limits_resolved.file_interpro_xml_bytes_max,
                label="InterPro XML file",
            )
        return cls(
            snapshot=_InterProSnapshot(
                file_protein2ipr=file_protein2ipr,
                file_interpro_xml=file_interpro_xml,
            ),
            limits=limits_resolved,
        )

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        kind_input_id: InterProInputIdKind,
    ) -> InterProSelection:
        """Create a selection for one normalized UniProt ID set.

        Args:
            ids: UniProt accessions. Empty values are discarded and duplicate
                normalized IDs collapse to one input row.
            kind_input_id: Input namespace. The supported value is `"uniprot"`.

        Returns:
            A selection handle whose mapping output includes `InputId` and
            `KindInputId` provenance columns.

        Raises:
            ValueError: If the namespace is unsupported or the normalized ID
                count exceeds the configured limit.

        Examples:
            Normalize a UniProt entry label before selecting its mapping rows:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(
            ...     ["sp|P12345|TEST_HUMAN"], kind_input_id="uniprot"
            ... )
            >>> selection.extract_mapping().select(
            ...     "InputId", "InterProId", "MemberDbId"
            ... ).rows()
            [('P12345', 'IPR000001', 'PF00051'), ('P12345', 'IPR000001', 'SM00130')]
        """
        validate_kind_input_id(kind_input_id)
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return InterProSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            kind_input_id=kind_input_id,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
        *,
        kind_input_id: InterProInputIdKind,
    ) -> InterProSelection:
        """Create a selection that preserves caller-defined groups.

        Args:
            group_to_ids: Mapping from unique, non-empty group IDs to UniProt
                accessions.
            kind_input_id: Input namespace. The supported value is `"uniprot"`.

        Returns:
            A selection handle whose mapping and unmapped outputs retain
            `GroupId`.

        Raises:
            ValueError: If group IDs or the namespace are invalid, or a group
                or input-count limit is exceeded.

        Examples:
            Preserve comparison labels on every selected mapping row:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"], "down": ["Q9Y243"]},
            ...     kind_input_id="uniprot",
            ... )
            >>> selection.extract_mapping().select(
            ...     "GroupId", "UniProtId"
            ... ).unique().sort("GroupId").rows()
            [('down', 'Q9Y243'), ('up', 'P12345')]
        """
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
        return InterProSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            kind_input_id=kind_input_id,
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full row-level UniProt-to-InterPro mapping.

        Returns:
            A cached frame that preserves member-database rows and positional
            coordinates. `InterProType` and `MemberDb` remain null when no XML
            source was configured.

        Examples:
            Retain each member-database match and its coordinates:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> db.extract_mapping().select(
            ...     "UniProtId", "MemberDbId", "Start", "End"
            ... ).head(2).rows()
            [('P12345', 'PF00051', 10, 80), ('P12345', 'SM00130', 12, 76)]
        """
        if self._df_mapping is None:
            self._df_mapping = read_mapping_frame(
                self.snapshot.file_protein2ipr,
                df_interpro_entry=self.xml_frame("entry"),
                df_interpro_member=self.xml_frame("member"),
            )
        return self._df_mapping

    def build_tidy(
        self,
        *,
        config: InterProTidyConfig = "mapping",
    ) -> InterProTidyDataset:
        """Build one configured lazy tidy dataset.

        Args:
            config: `"mapping"` builds the canonical `mapping.parquet` plan;
                `"pfam"` builds `protein_term`, `term`, and `term_xref` plans.

        Returns:
            A tidy dataset whose lazy frames and asset contract match `config`.

        Raises:
            ValueError: If `config` is unknown, Pfam XML is absent, or Pfam
                snapshot metadata and raw mappings fail strict validation.

        Notes:
            The mapping configuration accepts missing XML and leaves its
            enrichment fields null. The Pfam configuration requires both raw
            files under the same `<version>/raw/` directory and does not depend
            on a previously written canonical mapping.

        Examples:
            The compact ``108.0`` snapshot below must include an INTERPRO
            release record whose version is ``108.0``, matching its versioned
            parent directory.

            Build the default mapping plan:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     file_interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['mapping']

            Build the compact Pfam plan from the same raw snapshot:

            >>> sorted(db.build_tidy(config="pfam").frames)
            ['protein_term', 'term', 'term_xref']
        """
        self._validate_tidy_config(config)
        if config == "pfam":
            return build_pfam_tidy_dataset(
                file_protein2ipr=self.snapshot.file_protein2ipr,
                file_interpro_xml=self._require_interpro_xml(
                    action="build the Pfam tidy configuration"
                ),
            )

        return InterProTidyDataset(
            frames={
                "mapping": scan_mapping_frame(
                    self.snapshot.file_protein2ipr,
                    df_interpro_entry=self.xml_frame("entry"),
                    df_interpro_member=self.xml_frame("member"),
                )
            },
            source=self._tidy_sources(),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"interpro-mapping-{self.snapshot.file_protein2ipr.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
        )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        config: InterProTidyConfig = "mapping",
        should_write_manifest: bool = False,
        should_hash_sources: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write one configured InterPro tidy dataset.

        Args:
            dir_out: Destination directory for the configured parquet assets
                and optional manifest.
            config: `"mapping"` writes `mapping.parquet`; `"pfam"` writes the
                compact protein-term, term, and term-xref assets.
            should_write_manifest: Whether to write `manifest.json` and return
                its content in the report.
            should_hash_sources: Whether a requested manifest should contain
                SHA-256 values for every source file.
            should_hash_assets: Whether a requested manifest should contain
                SHA-256 values for every written parquet asset.

        Returns:
            A report describing the configured assets and optional manifest.

        Raises:
            ValueError: If `config` is unknown, Pfam XML is absent, or Pfam
                validation fails.

        Notes:
            Hash flags are publication costs whose values are observable only
            when `should_write_manifest=True`.

        Examples:
            The compact ``108.0`` snapshot below must include an INTERPRO
            release record whose version is ``108.0``, matching its versioned
            parent directory.

            Write the canonical mapping asset:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     file_interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> report = db.write_tidy("out/interpro")
            >>> [asset.path for asset in report.assets]
            ['mapping.parquet']

            Write the compact Pfam assets:

            >>> report = db.write_tidy("out/interpro-pfam", config="pfam")
            >>> [asset.path for asset in report.assets]
            ['protein_term.parquet', 'term.parquet', 'term_xref.parquet']
        """
        self._validate_tidy_config(config)
        if config == "pfam":
            dataset = build_pfam_tidy_dataset(
                file_protein2ipr=self.snapshot.file_protein2ipr,
                file_interpro_xml=self._require_interpro_xml(
                    action="write the Pfam tidy configuration"
                ),
                should_hash_sources=should_write_manifest and should_hash_sources,
            )
            return dataset.write(
                Path(dir_out),
                should_write_manifest=should_write_manifest,
                should_hash_assets=should_hash_assets,
            )

        dir_out = Path(dir_out)
        dir_out.mkdir(parents=True, exist_ok=True)
        file_out = dir_out / "mapping.parquet"
        scan_mapping_frame(
            self.snapshot.file_protein2ipr,
            df_interpro_entry=self.xml_frame("entry"),
            df_interpro_member=self.xml_frame("member"),
        ).sink_parquet(file_out)

        asset = TidyReportAsset(path="mapping.parquet", kind="canonical")
        asset_manifest = TidyManifestAsset(
            path=asset.path,
            kind=asset.kind,
            is_optional=asset.is_optional,
            sha256=calculate_file_sha256(file_out) if should_hash_assets else None,
        )
        manifest = (
            self._build_manifest(
                asset_manifest,
                should_hash_sources=should_hash_sources,
            )
            if should_write_manifest
            else None
        )
        if manifest is not None:
            (dir_out / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return TidyWriteReport(dir_out=dir_out, assets=(asset,), manifest=manifest)

    def _build_manifest(
        self,
        asset: TidyManifestAsset,
        *,
        should_hash_sources: bool,
    ) -> TidyManifest:
        timestamp = datetime.now(UTC)
        return {
            "build_id": f"interpro-mapping-{timestamp.strftime('%Y%m%dT%H%M%SZ')}",
            "schema_version": SCHEMA_VERSION,
            "generated_at": timestamp.isoformat().replace("+00:00", "Z"),
            "sources": [
                {
                    "path": source.path.as_posix(),
                    "bytes": source.path.stat().st_size,
                    "media_type": source.media_type,
                    **({"sha256": source.sha256} if source.sha256 is not None else {}),
                }
                for source in self._tidy_sources(
                    should_hash_sources=should_hash_sources
                )
            ],
            "assets": [asdict(asset)],
        }

    def xml_frame(self, frame_name: str) -> pl.DataFrame:
        """Read and cache one InterPro XML lookup frame.

        Args:
            frame_name: `"entry"` for entry types or `"member"` for
                InterPro-to-member-database relationships.

        Returns:
            The requested lookup frame. If no XML source was configured, the
            frame is empty but retains its expected schema.

        Raises:
            KeyError: If `frame_name` is not `"entry"` or `"member"`.

        Examples:
            Read the entry-type lookup from the configured XML:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     file_interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> db.xml_frame("entry").rows()
            [('IPR000001', 'Domain'), ('IPR000002', 'Homologous_superfamily')]
        """
        if self._frames_xml is None:
            self._frames_xml = read_interpro_xml_frames(self.snapshot.file_interpro_xml)
        return self._frames_xml[frame_name]

    def _tidy_sources(
        self,
        *,
        should_hash_sources: bool = False,
    ) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                path=self.snapshot.file_protein2ipr,
                media_type=MEDIA_TYPE_TSV_GZIP,
                sha256=calculate_file_sha256(self.snapshot.file_protein2ipr)
                if should_hash_sources
                else None,
            )
        ]
        if self.snapshot.file_interpro_xml is not None:
            sources.append(
                TidySource(
                    path=self.snapshot.file_interpro_xml,
                    media_type=MEDIA_TYPE_XML_GZIP,
                    sha256=calculate_file_sha256(self.snapshot.file_interpro_xml)
                    if should_hash_sources
                    else None,
                )
            )
        return tuple(sources)

    @staticmethod
    def _validate_tidy_config(config: str) -> None:
        if config not in {"mapping", "pfam"}:
            raise ValueError(
                f"InterPro tidy config must be one of ['mapping', 'pfam']: {config!r}"
            )

    def _require_interpro_xml(self, *, action: str) -> Path:
        if self.snapshot.file_interpro_xml is None:
            raise ValueError(f"InterPro XML file is required to {action}")
        return self.snapshot.file_interpro_xml


@dataclass(slots=True)
class InterProSelection:
    """Materialize one single or grouped InterPro mapping query.

    Instances are created by `InterProDb.select_ids()` or
    `InterProDb.select_groups()`. Mapping and unmapped frames are cached
    independently after first extraction.

    Examples:
        Inspect the InterPro IDs retained by a selection:

        >>> db = InterProDb.from_mapping_files(
        ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
        ... )
        >>> selection = db.select_ids(["P12345"], kind_input_id="uniprot")
        >>> selection.extract_mapping().get_column("InterProId").unique().to_list()
        ['IPR000001']
    """

    dataset: InterProDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    kind_input_id: InterProInputIdKind
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> db.select_groups(
            ...     {"up": ["P12345"]}, kind_input_id="uniprot"
            ... ).is_grouped
            True
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        return ("GroupId",) if self.is_grouped else ()

    def extract_mapping(self) -> pl.DataFrame:
        """Extract mapping rows for the normalized selection.

        Returns:
            A frame prefixed by `InputId` and `KindInputId`; grouped selections
            additionally begin with `GroupId`. Member and positional rows stay
            expanded.

        Examples:
            Extract the normalized input alongside its member-database rows:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(["P12345"], kind_input_id="uniprot")
            >>> selection.extract_mapping().select(
            ...     "InputId", "MemberDbId"
            ... ).rows()
            [('P12345', 'PF00051'), ('P12345', 'SM00130')]
        """
        if self._df_mapping is None:
            self._df_mapping = select_mapping_frame(
                self.dataset.snapshot.file_protein2ipr,
                self._df_input_ids,
                kind_input_id=self.kind_input_id,
                cols_group_id=self._col_group_id,
                df_interpro_entry=self.dataset.xml_frame("entry"),
                df_interpro_member=self.dataset.xml_frame("member"),
            )
        return self._df_mapping

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no InterPro mapping row.

        Returns:
            `InputId` for a single selection, or `GroupId, InputId` for a
            grouped selection.

        Examples:
            Report an accession absent from the local snapshot:

            >>> db = InterProDb.from_mapping_files(
            ...     file_protein2ipr="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(["MISSING"], kind_input_id="uniprot")
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
