from __future__ import annotations

import json
import os
import shutil
import tempfile
from collections.abc import Iterable
from dataclasses import asdict
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Literal, cast

import polars as pl

from bioextract._shared import validate_file_size
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

from .constant import (
    COLS_IDMAPPING_SELECTED,
    MEDIA_TYPE_FLAT_FILE,
    MEDIA_TYPE_FLAT_FILE_GZIP,
    MEDIA_TYPE_PARQUET,
    MEDIA_TYPE_PARQUET_DATASET,
    MEDIA_TYPE_TSV,
    MEDIA_TYPE_TSV_GZIP,
    SCHEMA_VERSION_EGGNOG_XREF,
    SCHEMA_VERSION_SUBCELLULAR_LOCATION,
    SCHEMA_VERSION,
)
from .util import (
    filter_taxids,
    has_hive_parquet_candidates,
    normalize_taxids,
    read_eggnog_xref_frame,
    read_subcellular_location_frame,
    scan_eggnog_xref_tsv,
    scan_hive_mapping_dataset,
    scan_raw_idmapping_selected,
    scan_subcellular_location_tsv,
    validate_mapping_schema,
    write_eggnog_xref_tsv,
    write_subcellular_location_tsv,
)

__all__ = [
    "UniprotDb",
    "UniprotResourceLimits",
]


class UniprotTidyManifest(TidyManifest):
    """Manifest schema for a taxid-scoped UniProt mapping publication."""

    taxids: list[str]


class _UniprotMappingKind(StrEnum):
    RAW_TSV = "raw_tsv"
    RAW_TSV_GZIP = "raw_tsv_gzip"
    PARQUET = "parquet"
    HIVE_PARQUET = "hive_parquet"
    DAT = "dat"
    DAT_GZIP = "dat_gzip"


@dataclass(frozen=True, slots=True)
class UniprotResourceLimits:
    """Optional guards for UniProt source-file sizes.

    Attributes:
        file_idmapping_selected_bytes_max: Maximum size of a raw TSV or single
            parquet mapping file, in bytes. `None` disables the limit. Dataset
            directories are checked structurally instead.
        file_dat_bytes_max: Maximum size of a UniProtKB `.dat` or `.dat.gz`
            source, in bytes. `None` disables the limit.

    Examples:
        Reject an oversized flat-file snapshot before parsing it:

        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_dat = Path(dir_tmp) / "uniprot_sprot.dat"
        ...     _ = file_dat.write_text("ID   P53_HUMAN\\n", encoding="utf-8")
        ...     limits = UniprotResourceLimits(file_dat_bytes_max=1)
        ...     try:
        ...         UniprotDb.from_dat(
        ...             file_dat=file_dat, source_db="sprot", limits=limits
        ...         )
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_idmapping_selected_bytes_max: int | None = None
    file_dat_bytes_max: int | None = None


@dataclass(frozen=True, slots=True)
class _UniprotSnapshot:
    kind: _UniprotMappingKind
    file_idmapping_selected: Path | None = None
    file_dat: Path | None = None
    taxids: tuple[str, ...] = ()
    source_db: str | None = None


@dataclass(slots=True)
class UniprotDb:
    """Access one UniProt mapping or knowledge-base flat-file snapshot.

    `from_files()` creates an idmapping handle over raw TSV, one normalized
    parquet, or a hive-partitioned parquet dataset. `from_dat()` creates a
    mutually exclusive flat-file handle for curated eggNOG cross-references and
    subcellular-location comments. Calling a method from the other mode raises
    `ValueError` instead of silently interpreting the wrong resource.

    Construction is deliberately lightweight: paths and configured file sizes
    are validated, but mapping data and schemas are not scanned until requested.

    Examples:
        Extract the human accessions from an idmapping resource:

        >>> db = UniprotDb.from_files(
        ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
        ... )
        >>> db.with_taxids("9606").extract_mapping().select(
        ...     "UniProtId", "GeneId"
        ... ).rows()
        [('P04637', '7157'), ('Q9Y243', '10000')]

        Resolve a secondary accession to its eggNOG cross-references:

        >>> kb = UniprotDb.from_dat(
        ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
        ...     source_db="sprot",
        ... )
        >>> kb.select_eggnog_xref_ids(["Q11111"]).select(
        ...     "PrimaryUniProtId", "EggnogOgId"
        ... ).rows()
        [('P12345', 'ENOG502ABC'), ('P12345', 'KOG0001')]
    """

    snapshot: _UniprotSnapshot
    limits: UniprotResourceLimits = field(default_factory=UniprotResourceLimits)

    DEFAULT_RESOURCE_LIMITS = UniprotResourceLimits()

    @classmethod
    def from_files(
        cls,
        *,
        file_idmapping_selected: os.PathLike[str] | str,
        limits: UniprotResourceLimits | None = None,
    ) -> UniprotDb:
        """Create a dataset handle from raw or tidy UniProt mapping data.

        Args:
            file_idmapping_selected: Path to `idmapping_selected.tab(.gz)`, a
                normalized parquet file, or a hive parquet dataset directory.
            limits: Dataset-level resource limits. Size limits apply to file
                inputs; hive dataset directories are checked structurally only.

        Returns:
            A dataset handle that can be taxid-scoped and materialized later.

        Raises:
            FileNotFoundError: If the path does not exist.
            ValueError: If the path type is unsupported, a directory contains
                no parquet files, or a configured file-size limit is exceeded.

        Examples:
            Read mouse records from a compressed idmapping source:

            >>> db = UniprotDb.from_files(
            ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
            ... )
            >>> db.with_taxids("10090").extract_mapping().select(
            ...     "UniProtId", "TaxId"
            ... ).rows()
            [('P31750', '10090')]
        """
        path = Path(file_idmapping_selected)
        if not path.exists():
            raise FileNotFoundError(
                f"UniProt idmapping selected path not found: {path}"
            )

        limits_resolved = UniprotResourceLimits() if limits is None else limits
        kind = _infer_mapping_kind(path)
        if path.is_file():
            validate_file_size(
                file_path=path,
                size_max=limits_resolved.file_idmapping_selected_bytes_max,
                label="UniProt idmapping selected file",
            )

        return cls(
            snapshot=_UniprotSnapshot(file_idmapping_selected=path, kind=kind),
            limits=limits_resolved,
        )

    @classmethod
    def from_dat(
        cls,
        *,
        file_dat: os.PathLike[str] | str,
        source_db: str,
        limits: UniprotResourceLimits | None = None,
    ) -> UniprotDb:
        """Create a dataset handle from a UniProtKB flat file.

        Args:
            file_dat: Path to a UniProtKB `.dat` or `.dat.gz` flat file.
            source_db: Non-empty provenance label copied to extracted rows and
                tidy build IDs.
            limits: Optional flat-file size guard.

        Returns:
            A flat-file handle for eggNOG xref and subcellular-location APIs.

        Raises:
            FileNotFoundError: If `file_dat` does not exist.
            ValueError: If `source_db` is empty, the suffix is unsupported, or
                the file exceeds its configured size limit.

        Notes:
            This handle cannot serve normalized idmapping extraction or
            `write_tidy()`; use `from_files()` for those operations.

        Examples:
            Extract eggNOG identifiers from a UniProtKB flat file:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_db="sprot",
            ... )
            >>> db.extract_eggnog_xref().select(
            ...     "EggnogOgId", "EggnogLevel"
            ... ).head(2).rows()
            [('ENOG502ABC', 'Metazoa'), ('KOG0001', 'Eukaryota')]
        """
        source_db = str(source_db).strip()
        if not source_db:
            raise ValueError("UniProt source_db must be non-empty after normalization")

        path = Path(file_dat)
        if not path.exists():
            raise FileNotFoundError(f"UniProt flat file not found: {path}")
        limits_resolved = UniprotResourceLimits() if limits is None else limits
        kind = _infer_dat_kind(path)
        validate_file_size(
            file_path=path,
            size_max=limits_resolved.file_dat_bytes_max,
            label="UniProt flat file",
        )
        return cls(
            snapshot=_UniprotSnapshot(file_dat=path, kind=kind, source_db=source_db),
            limits=limits_resolved,
        )

    def with_taxids(self, *taxids: str | int) -> UniprotDb:
        """Create a taxid-scoped view of an idmapping resource.

        Args:
            *taxids: NCBI taxonomy IDs. Values are normalized to non-empty
                strings and must remain distinct after normalization.

        Returns:
            A new handle sharing the same source and limits with the requested
            taxid scope.

        Raises:
            ValueError: If this is a flat-file handle, a taxid normalizes to an
                empty value, or normalized taxids are duplicated.

        Examples:
            Limit extraction to human records:

            >>> db = UniprotDb.from_files(
            ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
            ... )
            >>> db.with_taxids("9606").extract_mapping().select(
            ...     "UniProtId", "TaxId"
            ... ).rows()
            [('P04637', '9606'), ('Q9Y243', '9606')]
        """
        self._require_idmapping_snapshot("scope UniProt idmapping by taxid")
        return UniprotDb(
            snapshot=_UniprotSnapshot(
                file_idmapping_selected=self.snapshot.file_idmapping_selected,
                kind=self.snapshot.kind,
                taxids=normalize_taxids(taxids),
            ),
            limits=self.limits,
        )

    def validate_schema(self) -> None:
        """Validate the normalized idmapping schema without collecting rows.

        Raises:
            ValueError: If this is a flat-file handle or required mapping
                columns are missing.

        Examples:
            Normalize a compact raw fixture, then validate the resulting
            parquet without collecting it:

            >>> from pathlib import Path
            >>> from tempfile import TemporaryDirectory
            >>> raw_db = UniprotDb.from_files(
            ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
            ... ).with_taxids("9606")
            >>> with TemporaryDirectory() as dir_out:
            ...     dir_normalized = Path(dir_out) / "normalized"
            ...     report = raw_db.write_tidy(dir_normalized)
            ...     file_mapping = dir_normalized / report.assets[0].path
            ...     db = UniprotDb.from_files(file_idmapping_selected=file_mapping)
            ...     db.validate_schema() is None
            True
        """
        self._require_idmapping_snapshot("validate UniProt idmapping schema")
        validate_mapping_schema(self._scan_mapping())

    def extract_mapping(self) -> pl.DataFrame:
        """Extract normalized idmapping rows for the current taxid scope.

        Returns:
            All canonical mapping columns in contract order. An unscoped handle
            returns all taxids.

        Raises:
            ValueError: If this is a flat-file handle or the mapping schema is
                incomplete.

        Examples:
            Extract the identifiers retained by a human taxid scope:

            >>> db = UniprotDb.from_files(
            ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
            ... )
            >>> db.with_taxids("9606").extract_mapping().select(
            ...     "UniProtId", "GeneId", "TaxId"
            ... ).rows()
            [('P04637', '7157', '9606'), ('Q9Y243', '10000', '9606')]
        """
        self._require_idmapping_snapshot("extract UniProt idmapping")
        lf_mapping = self._scan_mapping()
        validate_mapping_schema(lf_mapping)
        return (
            filter_taxids(lf_mapping, self.snapshot.taxids)
            .select(COLS_IDMAPPING_SELECTED)
            .collect()
        )

    def extract_eggnog_xref(self) -> pl.DataFrame:
        """Extract all eggNOG cross-reference rows from a flat-file handle.

        Returns:
            A canonical frame that retains primary and secondary UniProt
            accessions, eggNOG identifiers, and source provenance.

        Raises:
            ValueError: If this is an idmapping handle.

        Examples:
            Retain each eggNOG identifier and its taxonomic level:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_db="sprot",
            ... )
            >>> db.extract_eggnog_xref().select(
            ...     "UniProtId", "EggnogOgId", "EggnogLevel"
            ... ).head(2).rows()
            [('P12345', 'ENOG502ABC', 'Metazoa'), ('P12345', 'KOG0001', 'Eukaryota')]
        """
        self._require_dat_snapshot("extract UniProt eggNOG xrefs")
        return read_eggnog_xref_frame(
            self._required_path(self.snapshot.file_dat),
            source_db=self.snapshot.source_db or "",
        )

    def select_eggnog_xref_ids(self, ids: Iterable[str]) -> pl.DataFrame:
        """Extract eggNOG xrefs for selected UniProt accessions.

        Args:
            ids: UniProt accessions. Empty values are discarded and duplicates
                collapse before the flat file is scanned.

        Returns:
            The canonical eggNOG xref frame restricted to matching accessions.

        Raises:
            ValueError: If this is an idmapping handle.

        Examples:
            Restrict the flat-file scan to one secondary accession:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_db="sprot",
            ... )
            >>> db.select_eggnog_xref_ids(["Q11111"]).select(
            ...     "UniProtId", "PrimaryUniProtId", "EggnogOgId"
            ... ).rows()
            [('Q11111', 'P12345', 'ENOG502ABC'), ('Q11111', 'P12345', 'KOG0001')]
        """
        self._require_dat_snapshot("select UniProt eggNOG xrefs")
        input_ids = {str(input_id).strip() for input_id in ids if str(input_id).strip()}
        return read_eggnog_xref_frame(
            self._required_path(self.snapshot.file_dat),
            source_db=self.snapshot.source_db or "",
            input_ids=input_ids,
        )

    def extract_subcellular_location(self) -> pl.DataFrame:
        """Extract curated UniProtKB subcellular-location comments.

        Returns:
            One row per accession, location text, and evidence reference,
            preserving comment notes and source provenance. Missing comments
            do not produce negative-location rows.

        Raises:
            ValueError: If this is an idmapping handle.

        Examples:
            Read curated locations together with their ECO evidence:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot_subcellular.dat.gz",
            ...     source_db="sprot",
            ... )
            >>> db.extract_subcellular_location().select(
            ...     "UniProtId", "SubcellularLocation", "EvidenceCode"
            ... ).head(2).rows()
            [('P12345', 'Cytoplasm', 'ECO:0000269'), ('P12345', 'Nucleus', 'ECO:0000305')]
        """
        self._require_dat_snapshot("extract UniProt subcellular locations")
        return read_subcellular_location_frame(
            self._required_path(self.snapshot.file_dat),
            source_db=self.snapshot.source_db or "",
        )

    def write_eggnog_xref_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write flat-file eggNOG xrefs as a canonical parquet mapping.

        Args:
            dir_out: Destination directory for `mapping.parquet` and the
                optional manifest.
            should_write_manifest: Whether to write `manifest.json` and return
                its content in the report.
            should_hash_assets: Whether a requested manifest should contain the
                mapping SHA-256 value.

        Returns:
            A report describing the mapping asset and optional manifest.

        Raises:
            ValueError: If this is an idmapping handle.

        Examples:
            Write the xref mapping and inspect the published asset name:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_db="Swiss-Prot",
            ... )
            >>> report = db.write_eggnog_xref_tidy("out/uniprot-eggnog")
            >>> [asset.path for asset in report.assets]
            ['mapping.parquet']
        """
        self._require_dat_snapshot("write UniProt eggNOG xref tidy")
        file_dat = self._required_path(self.snapshot.file_dat)
        media_type = (
            MEDIA_TYPE_FLAT_FILE_GZIP
            if self.snapshot.kind == _UniprotMappingKind.DAT_GZIP
            else MEDIA_TYPE_FLAT_FILE
        )
        with tempfile.TemporaryDirectory(prefix="bioextract-uniprot-xref-") as dir_tmp:
            file_xref_tsv = Path(dir_tmp) / "mapping.tsv"
            write_eggnog_xref_tsv(
                file_dat,
                file_xref_tsv,
                source_db=self.snapshot.source_db or "",
            )
            dataset = TidyDataset(
                frames={"mapping": scan_eggnog_xref_tsv(file_xref_tsv)},
                source=TidySource(path=file_dat, media_type=media_type),
                schema_version=SCHEMA_VERSION_EGGNOG_XREF,
                build_id_prefix=f"uniprot-eggnog-xref-{self.snapshot.source_db}",
                assets=(
                    TidyAsset(
                        path="mapping.parquet", kind="canonical", frame_name="mapping"
                    ),
                ),
            )
            return dataset.write(
                Path(dir_out),
                should_write_manifest=should_write_manifest,
                should_hash_assets=should_hash_assets,
            )

    def write_subcellular_location_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write curated subcellular locations as canonical parquet.

        Args:
            dir_out: Destination directory for `data.parquet` and the optional
                manifest.
            should_write_manifest: Whether to write `manifest.json` and return
                its content in the report.
            should_hash_assets: Whether a requested manifest should contain the
                data-asset SHA-256 value.

        Returns:
            A report describing the data asset and optional manifest.

        Raises:
            ValueError: If this is an idmapping handle.

        Examples:
            Write curated locations and inspect the published asset name:

            >>> db = UniprotDb.from_dat(
            ...     file_dat="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_db="Swiss-Prot",
            ... )
            >>> report = db.write_subcellular_location_tidy("out/uniprot-subcell")
            >>> [asset.path for asset in report.assets]
            ['data.parquet']
        """
        self._require_dat_snapshot("write UniProt subcellular location tidy")
        file_dat = self._required_path(self.snapshot.file_dat)
        media_type = (
            MEDIA_TYPE_FLAT_FILE_GZIP
            if self.snapshot.kind == _UniprotMappingKind.DAT_GZIP
            else MEDIA_TYPE_FLAT_FILE
        )
        with tempfile.TemporaryDirectory(
            prefix="bioextract-uniprot-subcell-"
        ) as dir_tmp:
            file_subcell_tsv = Path(dir_tmp) / "data.tsv"
            write_subcellular_location_tsv(
                file_dat,
                file_subcell_tsv,
                source_db=self.snapshot.source_db or "",
            )
            dataset = TidyDataset(
                frames={"data": scan_subcellular_location_tsv(file_subcell_tsv)},
                source=TidySource(path=file_dat, media_type=media_type),
                schema_version=SCHEMA_VERSION_SUBCELLULAR_LOCATION,
                build_id_prefix=f"uniprot-subcellular-location-{self.snapshot.source_db}",
                assets=(
                    TidyAsset(path="data.parquet", kind="canonical", frame_name="data"),
                ),
            )
            return dataset.write(
                Path(dir_out),
                should_write_manifest=should_write_manifest,
                should_hash_assets=should_hash_assets,
            )

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_allow_all: bool = False,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
        policy_existing: Literal["error", "overwrite", "skip"] = "error",
        dir_tmp: os.PathLike[str] | str | None = None,
        level_compression: int | None = None,
        should_monitor_resources: bool = True,
        size_rss_stop_gb: float | None = 96,
        size_threads_stop: int | None = 160,
        count_d_state_stop: int = 3,
    ) -> TidyWriteReport:
        """Write normalized UniProt mapping data as parquet.

        Args:
            dir_out: Output directory.
            should_allow_all: Required when no taxids are selected, because all
                taxa may be very large.
            should_write_manifest: Whether to write `manifest.json`.
            should_hash_assets: Whether to calculate asset checksums in the
                manifest.
            policy_existing: How to handle a non-empty output directory:
                `error`, `overwrite`, or `skip`.
            dir_tmp: Optional scratch directory for staging output before it is
                published to `dir_out`.
            level_compression: Optional zstd compression level for
                `mapping.parquet`.
            should_monitor_resources: Whether to stop when current-process
                resource limits are exceeded.
            size_rss_stop_gb: Stop threshold for current-process RSS.
            size_threads_stop: Stop threshold for current-process thread count.
            count_d_state_stop: Number of consecutive D-state samples that
                trigger a stop.

        Returns:
            A write report with asset paths and optional manifest content.

        Raises:
            FileExistsError: If `dir_out` is non-empty and
                `policy_existing="error"`.
            ValueError: If this is a flat-file handle, the policy is invalid,
                an all-taxid write was not acknowledged, the schema is invalid,
                or CephFS staging requirements are not met.
            RuntimeError: If enabled RSS, thread, or D-state safety thresholds
                are exceeded.

        Notes:
            The mapping is written to a staging directory and published only
            after the parquet and optional manifest are complete. Set a stop
            threshold to `None` to disable that individual resource check.

        Examples:
            Write a taxid-scoped mapping and inspect the published asset:

            >>> db = UniprotDb.from_files(
            ...     file_idmapping_selected="fixtures/uniprot/idmapping_selected.tab.gz"
            ... ).with_taxids("9606")
            >>> report = db.write_tidy("out/uniprot-human")
            >>> [asset.path for asset in report.assets]
            ['mapping.parquet']
        """
        self._require_idmapping_snapshot("write UniProt idmapping tidy")
        if policy_existing not in {"error", "overwrite", "skip"}:
            raise ValueError(
                "policy_existing must be one of: 'error', 'overwrite', 'skip'"
            )
        if not self.snapshot.taxids and not should_allow_all:
            raise ValueError(
                "Writing all UniProt taxids requires should_allow_all=True"
            )

        dir_out = Path(dir_out)
        if should_allow_all and _is_ceph_path(dir_out):
            if dir_tmp is None:
                raise ValueError(
                    "Writing all UniProt taxids to /cephfs_data requires an "
                    "explicit local dir_tmp"
                )
            if _is_ceph_path(Path(dir_tmp)):
                raise ValueError(
                    "UniProt dir_tmp must not be under /cephfs_data for "
                    "all-taxid writes"
                )

        if dir_out.exists() and any(dir_out.iterdir()):
            if policy_existing == "error":
                raise FileExistsError(
                    f"UniProt tidy output directory is not empty: {dir_out}"
                )
            if policy_existing == "skip":
                return _build_existing_tidy_report(
                    dir_out=dir_out,
                    should_write_manifest=should_write_manifest,
                )

        lf_mapping = filter_taxids(self._scan_mapping(), self.snapshot.taxids)
        validate_mapping_schema(lf_mapping)
        monitor = _ResourceMonitor(
            should_monitor_resources=should_monitor_resources,
            size_rss_stop_gb=size_rss_stop_gb,
            size_threads_stop=size_threads_stop,
            count_d_state_stop=count_d_state_stop,
        )

        dir_tmp_parent = None if dir_tmp is None else Path(dir_tmp)
        if dir_tmp_parent is not None:
            dir_tmp_parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(
            prefix="bioextract-uniprot-",
            dir=None if dir_tmp_parent is None else dir_tmp_parent,
        ) as dir_stage_raw:
            dir_stage = Path(dir_stage_raw) / "tidy"
            dir_stage.mkdir(parents=True, exist_ok=True)

            monitor.check()
            file_out = dir_stage / "mapping.parquet"
            lf_mapping.select(COLS_IDMAPPING_SELECTED).sink_parquet(
                file_out,
                compression="zstd",
                compression_level=level_compression,
            )
            monitor.check()
            assets = (TidyReportAsset(path="mapping.parquet", kind="canonical"),)
            assets_manifest: list[TidyManifestAsset] = [
                TidyManifestAsset(
                    path=asset.path,
                    kind=asset.kind,
                    is_optional=asset.is_optional,
                    sha256=calculate_file_sha256(file_out)
                    if should_hash_assets
                    else None,
                )
                for asset in assets
            ]

            manifest = (
                self._build_manifest(assets_manifest) if should_write_manifest else None
            )
            if manifest is not None:
                (dir_stage / "manifest.json").write_text(
                    json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )

            _publish_tidy_dir(
                dir_stage=dir_stage,
                dir_out=dir_out,
            )
        return TidyWriteReport(dir_out=dir_out, assets=assets, manifest=manifest)

    def _scan_mapping(self) -> pl.LazyFrame:
        self._require_idmapping_snapshot("scan UniProt idmapping")
        match self.snapshot.kind:
            case _UniprotMappingKind.RAW_TSV | _UniprotMappingKind.RAW_TSV_GZIP:
                return scan_raw_idmapping_selected(
                    self._required_path(self.snapshot.file_idmapping_selected)
                )
            case _UniprotMappingKind.PARQUET:
                return pl.scan_parquet(
                    self._required_path(self.snapshot.file_idmapping_selected)
                )
            case _UniprotMappingKind.HIVE_PARQUET:
                return scan_hive_mapping_dataset(
                    self._required_path(self.snapshot.file_idmapping_selected)
                )
            case _:
                raise ValueError(
                    "Cannot scan UniProt idmapping from flat-file snapshot"
                )

    def _build_manifest(
        self,
        assets: list[TidyManifestAsset],
    ) -> UniprotTidyManifest:
        timestamp = datetime.now(UTC)
        source: dict[str, str | int] = {
            "path": self._required_path(
                self.snapshot.file_idmapping_selected
            ).as_posix(),
            "media_type": _media_type_for_kind(self.snapshot.kind),
        }
        file_idmapping_selected = self._required_path(
            self.snapshot.file_idmapping_selected
        )
        if file_idmapping_selected.is_file():
            source["bytes"] = file_idmapping_selected.stat().st_size
        return {
            "build_id": f"uniprot-idmapping-selected-{timestamp.strftime('%Y%m%dT%H%M%SZ')}",
            "schema_version": SCHEMA_VERSION,
            "generated_at": timestamp.isoformat().replace("+00:00", "Z"),
            "taxids": list(self.snapshot.taxids),
            "sources": [source],
            "assets": [asdict(asset) for asset in assets],
        }

    def _require_idmapping_snapshot(self, action: str) -> None:
        if self.snapshot.kind in {
            _UniprotMappingKind.DAT,
            _UniprotMappingKind.DAT_GZIP,
        }:
            raise ValueError(f"Cannot {action} from a UniProt flat-file snapshot")
        if self.snapshot.file_idmapping_selected is None:
            raise ValueError(f"Cannot {action}: idmapping selected path is missing")

    def _require_dat_snapshot(self, action: str) -> None:
        if self.snapshot.kind not in {
            _UniprotMappingKind.DAT,
            _UniprotMappingKind.DAT_GZIP,
        }:
            raise ValueError(f"Cannot {action} from a UniProt idmapping snapshot")
        if self.snapshot.file_dat is None:
            raise ValueError(f"Cannot {action}: UniProt flat-file path is missing")

    @staticmethod
    def _required_path(path: Path | None) -> Path:
        if path is None:
            raise ValueError("Required UniProt resource path is missing")
        return path


def _infer_mapping_kind(path: Path) -> _UniprotMappingKind:
    if path.is_dir():
        if not has_hive_parquet_candidates(path):
            raise ValueError(
                f"UniProt hive parquet dataset contains no parquet files: {path}"
            )
        return _UniprotMappingKind.HIVE_PARQUET
    match path.name:
        case name if name.endswith(".tab.gz"):
            return _UniprotMappingKind.RAW_TSV_GZIP
        case name if name.endswith(".tab"):
            return _UniprotMappingKind.RAW_TSV
        case name if name.endswith(".parquet"):
            return _UniprotMappingKind.PARQUET
        case _:
            raise ValueError(
                f"Unsupported UniProt idmapping selected input type: {path}"
            )


def _infer_dat_kind(path: Path) -> _UniprotMappingKind:
    match path.name:
        case name if name.endswith(".dat.gz"):
            return _UniprotMappingKind.DAT_GZIP
        case name if name.endswith(".dat"):
            return _UniprotMappingKind.DAT
        case _:
            raise ValueError(f"Unsupported UniProt flat-file input type: {path}")


def _publish_tidy_dir(
    *,
    dir_stage: Path,
    dir_out: Path,
) -> None:
    dir_out.parent.mkdir(parents=True, exist_ok=True)
    if dir_out.exists():
        shutil.rmtree(dir_out)
    shutil.move(dir_stage.as_posix(), dir_out.as_posix())


def _build_existing_tidy_report(
    *,
    dir_out: Path,
    should_write_manifest: bool,
) -> TidyWriteReport:
    manifest = None
    file_manifest = dir_out / "manifest.json"
    if should_write_manifest and file_manifest.exists():
        manifest = cast(
            UniprotTidyManifest,
            json.loads(file_manifest.read_text(encoding="utf-8")),
        )
        assets = _extract_report_assets_from_manifest(manifest)
    else:
        assets: tuple[TidyReportAsset, ...] = (
            TidyReportAsset(path="mapping.parquet", kind="canonical"),
        )
    return TidyWriteReport(dir_out=dir_out, assets=assets, manifest=manifest)


def _extract_report_assets_from_manifest(
    manifest: TidyManifest,
) -> tuple[TidyReportAsset, ...]:
    return tuple(
        TidyReportAsset(
            path=str(asset["path"]),
            kind=str(asset["kind"]),
            is_optional=bool(asset["is_optional"]),
        )
        for asset in manifest["assets"]
    )


def _is_ceph_path(path: Path) -> bool:
    return path.as_posix() == "/cephfs_data" or path.as_posix().startswith(
        "/cephfs_data/"
    )


@dataclass(slots=True)
class _ResourceSample:
    size_rss_mb: float | None
    size_threads: int | None
    count_d_state_threads: int | None


@dataclass(slots=True)
class _ResourceMonitor:
    should_monitor_resources: bool
    size_rss_stop_gb: float | None
    size_threads_stop: int | None
    count_d_state_stop: int
    _count_d_state_consecutive: int = 0
    _size_threads_baseline: int | None = None

    def check(self) -> None:
        if not self.should_monitor_resources:
            return
        sample = _sample_process_resources()
        if (
            self.size_rss_stop_gb is not None
            and sample.size_rss_mb is not None
            and sample.size_rss_mb > self.size_rss_stop_gb * 1024
        ):
            raise RuntimeError(
                "UniProt tidy write exceeded RSS stop threshold: "
                f"{sample.size_rss_mb:.1f} MiB > {self.size_rss_stop_gb} GiB"
            )
        self._check_threads(sample.size_threads)
        if (
            sample.count_d_state_threads is not None
            and sample.count_d_state_threads > 0
        ):
            self._count_d_state_consecutive += 1
        else:
            self._count_d_state_consecutive = 0
        if self._count_d_state_consecutive >= self.count_d_state_stop:
            raise RuntimeError(
                "UniProt tidy write observed D-state threads in "
                f"{self._count_d_state_consecutive} consecutive samples"
            )

    def _check_threads(self, size_threads: int | None) -> None:
        if self.size_threads_stop is None or size_threads is None:
            return
        if self._size_threads_baseline is None:
            self._size_threads_baseline = size_threads
            if size_threads <= self.size_threads_stop or self.size_threads_stop < 64:
                self._raise_if_threads_exceeded(size_threads)
            return
        size_threads_stop_effective = self.size_threads_stop
        if self._size_threads_baseline > self.size_threads_stop:
            size_threads_stop_effective = self._size_threads_baseline + max(
                16,
                self.size_threads_stop // 4,
            )
        if size_threads > size_threads_stop_effective:
            raise RuntimeError(
                "UniProt tidy write exceeded thread stop threshold: "
                f"{size_threads} > {size_threads_stop_effective}"
            )

    def _raise_if_threads_exceeded(self, size_threads: int) -> None:
        if self.size_threads_stop is not None and size_threads > self.size_threads_stop:
            raise RuntimeError(
                "UniProt tidy write exceeded thread stop threshold: "
                f"{size_threads} > {self.size_threads_stop}"
            )


def _sample_process_resources() -> _ResourceSample:
    status = _read_proc_self_status()
    return _ResourceSample(
        size_rss_mb=_parse_status_size_mb(status.get("VmRSS")),
        size_threads=_parse_status_int(status.get("Threads")),
        count_d_state_threads=_count_proc_self_threads_by_state("D"),
    )


def _read_proc_self_status() -> dict[str, str]:
    try:
        lines = Path("/proc/self/status").read_text(encoding="utf-8").splitlines()
    except OSError:
        return {}
    status: dict[str, str] = {}
    for line in lines:
        key, sep, value = line.partition(":")
        if sep:
            status[key] = value.strip()
    return status


def _parse_status_size_mb(value: str | None) -> float | None:
    if value is None:
        return None
    parts = value.split()
    if not parts:
        return None
    try:
        size_kb = float(parts[0])
    except ValueError:
        return None
    return size_kb / 1024


def _parse_status_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value.split()[0])
    except (IndexError, ValueError):
        return None


def _count_proc_self_threads_by_state(state: str) -> int | None:
    dir_task = Path("/proc/self/task")
    try:
        files_stat = list(dir_task.glob("*/stat"))
    except OSError:
        return None
    count = 0
    for file_stat in files_stat:
        try:
            text = file_stat.read_text(encoding="utf-8")
        except OSError:
            continue
        state_thread = _parse_proc_stat_state(text)
        if state_thread == state:
            count += 1
    return count


def _parse_proc_stat_state(text: str) -> str | None:
    idx_close = text.rfind(")")
    if idx_close < 0:
        return None
    fields_after_name = text[idx_close + 1 :].strip().split()
    if not fields_after_name:
        return None
    return fields_after_name[0]


def _media_type_for_kind(kind: _UniprotMappingKind) -> str:
    match kind:
        case _UniprotMappingKind.RAW_TSV_GZIP:
            return MEDIA_TYPE_TSV_GZIP
        case _UniprotMappingKind.RAW_TSV:
            return MEDIA_TYPE_TSV
        case _UniprotMappingKind.PARQUET:
            return MEDIA_TYPE_PARQUET
        case _UniprotMappingKind.HIVE_PARQUET:
            return MEDIA_TYPE_PARQUET_DATASET
        case _UniprotMappingKind.DAT:
            return MEDIA_TYPE_FLAT_FILE
        case _UniprotMappingKind.DAT_GZIP:
            return MEDIA_TYPE_FLAT_FILE_GZIP
        case _:
            raise ValueError(f"Unsupported UniProt mapping kind: {kind!r}")
