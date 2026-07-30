from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import polars as pl

from bioextract._publication import ParquetWriteResult
from bioextract._tidy import (
    TidyAsset,
    TidyDataset,
    TidySource,
)

from .constant import (
    COLS_IDMAPPING_SELECTED,
    MEDIA_TYPE_FLAT_FILE,
    MEDIA_TYPE_FLAT_FILE_GZIP,
    MEDIA_TYPE_PARQUET,
    MEDIA_TYPE_PARQUET_DATASET,
    MEDIA_TYPE_TSV,
    MEDIA_TYPE_TSV_GZIP,
    SCHEMA_VERSION,
    SCHEMA_VERSION_EGGNOG_XREF,
    SCHEMA_VERSION_SUBCELLULAR_LOCATION,
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
    "UniProtDatabase",
]


class _UniprotMappingKind(StrEnum):
    RAW_TSV = "raw_tsv"
    RAW_TSV_GZIP = "raw_tsv_gzip"
    PARQUET = "parquet"
    HIVE_PARQUET = "hive_parquet"
    DAT = "dat"
    DAT_GZIP = "dat_gzip"


@dataclass(frozen=True, slots=True)
class _UniprotSnapshot:
    kind: _UniprotMappingKind
    file_idmapping_selected: Path | None = None
    file_dat: Path | None = None
    taxids: tuple[str, ...] = ()
    source_db: str | None = None


@dataclass(slots=True)
class UniProtDatabase:
    """Access one UniProt mapping or knowledge-base flat-file snapshot.

    `from_files()` creates an idmapping handle over raw TSV, one normalized
    parquet, or a hive-partitioned parquet dataset. `from_dat()` creates a
    mutually exclusive flat-file handle for curated eggNOG cross-references and
    subcellular-location comments. Calling a method from the other mode raises
    `ValueError` instead of silently interpreting the wrong resource.

    Construction is deliberately lightweight: paths are validated, but mapping
    data and schemas are not scanned until requested.

    Examples:
        Extract the human accessions from an idmapping resource:

        >>> db = UniProtDatabase.from_files(
        ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
        ... )
        >>> db.with_taxids("9606").extract_mapping().select(
        ...     "UniProtId", "GeneId"
        ... ).rows()
        [('P04637', '7157'), ('Q9Y243', '10000')]

        Resolve a secondary accession to its eggNOG cross-references:

        >>> kb = UniProtDatabase.from_dat(
        ...     path="fixtures/uniprot/uniprot_sprot.dat.gz",
        ...     source_database="sprot",
        ... )
        >>> kb.select_eggnog_xref_ids(["Q11111"]).select(
        ...     "PrimaryUniProtId", "EggnogOgId"
        ... ).rows()
        [('P12345', 'ENOG502ABC'), ('P12345', 'KOG0001')]
    """

    snapshot: _UniprotSnapshot

    @classmethod
    def from_files(
        cls,
        *,
        id_mapping: os.PathLike[str] | str,
    ) -> UniProtDatabase:
        """Create a dataset handle from raw or tidy UniProt mapping data.

        Args:
            id_mapping: Path to `idmapping_selected.tab(.gz)`, a
                normalized parquet file, or a hive parquet dataset directory.

        Returns:
            A dataset handle that can be taxid-scoped and materialized later.

        Raises:
            FileNotFoundError: If the path does not exist.
            ValueError: If the path type is unsupported, a directory contains
                no parquet files.

        Examples:
            Read mouse records from a compressed idmapping source:

            >>> db = UniProtDatabase.from_files(
            ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
            ... )
            >>> db.with_taxids("10090").extract_mapping().select(
            ...     "UniProtId", "TaxId"
            ... ).rows()
            [('P31750', '10090')]
        """
        path = Path(id_mapping)
        if not path.exists():
            raise FileNotFoundError(
                f"UniProt idmapping selected path not found: {path}"
            )

        kind = _infer_mapping_kind(path)
        return cls(snapshot=_UniprotSnapshot(file_idmapping_selected=path, kind=kind))

    @classmethod
    def from_dat(
        cls,
        path: os.PathLike[str] | str,
        *,
        source_database: str,
    ) -> UniProtDatabase:
        """Create a dataset handle from a UniProtKB flat file.

        Args:
            path: UniProtKB `.dat` or `.dat.gz` flat file.
            source_database: Non-empty provenance label copied to extracted rows and
                tidy build IDs.

        Returns:
            A flat-file handle for eggNOG xref and subcellular-location APIs.

        Raises:
            FileNotFoundError: If `path` does not exist.
            ValueError: If `source_database` is empty or the suffix is unsupported.

        Notes:
            This handle cannot serve normalized idmapping extraction; use
            `from_files()` for that operation.

        Examples:
            Extract eggNOG identifiers from a UniProtKB flat file:

            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_database="sprot",
            ... )
            >>> db.extract_eggnog_xref().select(
            ...     "EggnogOgId", "EggnogLevel"
            ... ).head(2).rows()
            [('ENOG502ABC', 'Metazoa'), ('KOG0001', 'Eukaryota')]
        """
        source_db = str(source_database).strip()
        if not source_db:
            raise ValueError(
                "UniProt source_database must be non-empty after normalization"
            )

        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"UniProt flat file not found: {path}")
        kind = _infer_dat_kind(path)
        return cls(
            snapshot=_UniprotSnapshot(file_dat=path, kind=kind, source_db=source_db)
        )

    def with_taxids(self, *taxids: str | int) -> UniProtDatabase:
        """Create a taxid-scoped view of an idmapping resource.

        Args:
            *taxids: NCBI taxonomy IDs. Values are normalized to non-empty
                strings and must remain distinct after normalization.

        Returns:
            A new handle sharing the same source with the requested taxid scope.

        Raises:
            ValueError: If this is a flat-file handle, a taxid normalizes to an
                empty value, or normalized taxids are duplicated.

        Examples:
            Limit extraction to human records:

            >>> db = UniProtDatabase.from_files(
            ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
            ... )
            >>> db.with_taxids("9606").extract_mapping().select(
            ...     "UniProtId", "TaxId"
            ... ).rows()
            [('P04637', '9606'), ('Q9Y243', '9606')]
        """
        self._require_idmapping_snapshot("scope UniProt idmapping by taxid")
        return UniProtDatabase(
            snapshot=_UniprotSnapshot(
                file_idmapping_selected=self.snapshot.file_idmapping_selected,
                kind=self.snapshot.kind,
                taxids=normalize_taxids(taxids),
            ),
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
            >>> raw_db = UniProtDatabase.from_files(
            ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
            ... ).with_taxids("9606")
            >>> with TemporaryDirectory() as dir_out:
            ...     file_mapping = Path(dir_out) / "normalized.parquet"
            ...     raw_db.write_parquet(file_mapping)
            ...     db = UniProtDatabase.from_files(id_mapping=file_mapping)
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

            >>> db = UniProtDatabase.from_files(
            ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
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

            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_database="sprot",
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

            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_database="sprot",
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

            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot_subcellular.dat.gz",
            ...     source_database="sprot",
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

    def write_eggnog_xref_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Stream UniProtKB eggNOG cross-references to one atomic Parquet.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot.dat.gz",
            ...     source_database="sprot",
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     db.write_eggnog_xref_parquet(
            ...         Path(dir_out) / "uniprot_eggnog.parquet"
            ...     ).schema_version
            'uniprot-eggnog-xref-v0.1'
        """
        self._require_dat_snapshot("write UniProt eggNOG xref parquet")
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
                resource_name="uniprot",
                frames={"mapping": scan_eggnog_xref_tsv(file_xref_tsv)},
                source=TidySource(path=file_dat, media_type=media_type),
                schema_version=SCHEMA_VERSION_EGGNOG_XREF,
                build_id_prefix=f"uniprot-eggnog-xref-{self.snapshot.source_db}",
                assets=(
                    TidyAsset(
                        path="mapping.parquet",
                        kind="canonical",
                        frame_name="mapping",
                    ),
                ),
            )
            return dataset.write_parquet(path, if_exists=if_exists)

    def write_subcellular_location_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Stream curated subcellular locations to one atomic Parquet.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = UniProtDatabase.from_dat(
            ...     path="fixtures/uniprot/uniprot_sprot_subcellular.dat.gz",
            ...     source_database="sprot",
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     db.write_subcellular_location_parquet(
            ...         Path(dir_out) / "uniprot_subcellular_location.parquet"
            ...     ).schema_version
            'uniprot-subcellular-location-v0.1'
        """
        self._require_dat_snapshot("write UniProt subcellular location parquet")
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
                resource_name="uniprot",
                frames={"data": scan_subcellular_location_tsv(file_subcell_tsv)},
                source=TidySource(path=file_dat, media_type=media_type),
                schema_version=SCHEMA_VERSION_SUBCELLULAR_LOCATION,
                build_id_prefix=(
                    f"uniprot-subcellular-location-{self.snapshot.source_db}"
                ),
                assets=(
                    TidyAsset(
                        path="data.parquet",
                        kind="canonical",
                        frame_name="data",
                    ),
                ),
            )
            return dataset.write_parquet(path, if_exists=if_exists)

    def write_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        allow_all_taxa: bool = False,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Stream the current idmapping scope to one atomic Parquet file.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = UniProtDatabase.from_files(
            ...     id_mapping="fixtures/uniprot/idmapping_selected.tab.gz"
            ... ).with_taxids("9606")
            >>> with TemporaryDirectory() as dir_out:
            ...     db.write_parquet(
            ...         Path(dir_out) / "uniprot.parquet"
            ...     ).schema_version
            'uniprot-idmapping-selected-v0.1'
        """
        self._require_idmapping_snapshot("write UniProt idmapping parquet")
        if not self.snapshot.taxids and not allow_all_taxa:
            raise ValueError("Writing all UniProt taxids requires allow_all_taxa=True")
        source_path = self._required_path(self.snapshot.file_idmapping_selected)
        frame = filter_taxids(self._scan_mapping(), self.snapshot.taxids).select(
            COLS_IDMAPPING_SELECTED
        )
        validate_mapping_schema(frame)
        dataset = TidyDataset(
            resource_name="uniprot",
            frames={"mapping": frame},
            source=TidySource(
                path=source_path,
                media_type=_media_type_for_kind(self.snapshot.kind),
            ),
            schema_version=SCHEMA_VERSION,
            build_id_prefix="uniprot-id-mapping",
            assets=(
                TidyAsset(
                    path="mapping.parquet",
                    kind="canonical",
                    frame_name="mapping",
                ),
            ),
        )
        return dataset.write_parquet(
            path,
            if_exists=if_exists,
            preserve_source_headers=True,
        )

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
