from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import duckdb
import polars as pl

from bioextract._publication import DuckDBWriteResult
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError, IntegrityError

from ._knowledgebase import write_knowledgebase
from ._query import UniProtSelection, make_selection, validate_publication
from .constant import (
    COLS_IDMAPPING_SELECTED,
    IDMAPPING_SOURCE_SCHEMA_PROFILE,
    MEDIA_TYPE_PARQUET,
    MEDIA_TYPE_PARQUET_DATASET,
    MEDIA_TYPE_TSV,
    MEDIA_TYPE_TSV_GZIP,
    SCHEMA_VERSION,
)
from .util import (
    filter_taxids,
    has_hive_parquet_candidates,
    normalize_taxids,
    scan_hive_mapping_dataset,
    scan_raw_idmapping_selected,
    validate_mapping_schema,
)

__all__ = ["UniProtDatabase"]


@dataclass(slots=True)
class UniProtDatabase:
    """Access one UniProt idmapping product or UniProtKB publication.

    Examples:
        >>> UniProtDatabase.from_idmapping  # doctest: +ELLIPSIS
        <bound method ...>
    """

    _mode: str
    _source_path: Path | None = None
    _mapping_kind: str | None = None
    _taxon_ids: tuple[str, ...] = ()
    _entries: Path | None = None
    _canonical_sequences: Path | None = None
    _isoform_sequences: Path | None = None
    _release_version: str | None = None
    _duckdb_path: Path | None = None
    _duckdb_identity: tuple[int, int, int, int, int] | None = None

    @classmethod
    def from_idmapping(
        cls,
        path: os.PathLike[str] | str,
        *,
        release_version: str | None = None,
    ) -> UniProtDatabase:
        """Open the separate lazy idmapping product.

        Examples:
            >>> UniProtDatabase.from_idmapping("mapping.parquet")  # doctest: +SKIP
            UniProtDatabase(...)
        """
        source = Path(path).resolve()
        if not source.exists():
            raise FileNotFoundError(f"UniProt idmapping path not found: {source}")
        kind = _mapping_kind(source)
        return cls(
            "idmapping",
            _source_path=source,
            _mapping_kind=kind,
            _release_version=_normalize_version(release_version),
        )

    @classmethod
    def from_knowledgebase(
        cls,
        *,
        entries: os.PathLike[str] | str,
        canonical_sequences: os.PathLike[str] | str | None = None,
        isoform_sequences: os.PathLike[str] | str | None = None,
        release_version: str | None = None,
    ) -> UniProtDatabase:
        """Declare exact UniProtKB input roles.

        Examples:
            >>> UniProtDatabase.from_knowledgebase(entries="uniprot.dat.gz")  # doctest: +SKIP
            UniProtDatabase(...)
        """
        entry_path = _required_file(entries, "entries")
        canonical_path = (
            None
            if canonical_sequences is None
            else _required_file(canonical_sequences, "canonical_sequences")
        )
        isoform_path = (
            None
            if isoform_sequences is None
            else _required_file(isoform_sequences, "isoform_sequences")
        )
        return cls(
            "knowledgebase_source",
            _entries=entry_path,
            _canonical_sequences=canonical_path,
            _isoform_sequences=isoform_path,
            _release_version=_normalize_version(release_version),
        )

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> UniProtDatabase:
        """Open and validate a UniProtKB publication.

        Examples:
            >>> UniProtDatabase.from_duckdb("uniprot.duckdb")  # doctest: +SKIP
            UniProtDatabase(...)
        """
        publication = Path(path).resolve()
        if not publication.is_file():
            raise FileNotFoundError(
                f"UniProt DuckDB publication not found: {publication}"
            )
        identity_before = _file_identity(publication)
        profile, metadata = validate_publication(publication)
        identity_after = _file_identity(publication)
        if identity_after != identity_before:
            raise IntegrityError("UniProt publication changed during validation")
        return cls(
            f"{profile}_publication",
            _duckdb_path=publication,
            _duckdb_identity=identity_after,
            _release_version=metadata.get("bioextract.release_version"),
        )

    @property
    def release_version(self) -> str | None:
        """Return the caller-supplied official release identity.

        Examples:
            >>> database.release_version  # doctest: +SKIP
            '2026_01'
        """
        return self._release_version

    def scan_mapping(
        self, *, taxon_ids: Iterable[str | int] | None = None
    ) -> pl.LazyFrame:
        """Return a lazy idmapping scan.

        Examples:
            >>> database.scan_mapping(taxon_ids=["9606"])  # doctest: +SKIP
            <LazyFrame ...>
        """
        if self._mode not in {"idmapping", "idmapping_publication"}:
            raise CapabilityError("UniProt publication has no idmapping capability")
        normalized = normalize_taxids(tuple(taxon_ids or ()))
        frame = self._scan_mapping(normalized)
        validate_mapping_schema(frame)
        if self._mode == "idmapping_publication":
            return frame.select(COLS_IDMAPPING_SELECTED)
        return filter_taxids(frame, normalized).select(COLS_IDMAPPING_SELECTED)

    def read_mapping(
        self,
        *,
        taxon_ids: Iterable[str | int] | None = None,
        allow_all_taxa: bool = False,
    ) -> pl.DataFrame:
        """Materialize an explicitly scoped idmapping read.

        Examples:
            >>> database.read_mapping(taxon_ids=["9606"])  # doctest: +SKIP
            shape: (...)
        """
        normalized = normalize_taxids(tuple(taxon_ids or ()))
        if not normalized and not allow_all_taxa:
            raise ValueError(
                "Unscoped eager idmapping reads require allow_all_taxa=True"
            )
        return self.scan_mapping(taxon_ids=normalized).collect()

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        taxon_ids: Iterable[str | int] | None = None,
        allow_all_taxa: bool = False,
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Publish idmapping or the declared UniProtKB roles as atomic DuckDB.

        Examples:
            >>> database.write_duckdb("mapping.duckdb", taxon_ids=["9606"])  # doctest: +SKIP
            DuckDBWriteResult(...)
        """
        if self._mode == "knowledgebase_source":
            if taxon_ids is not None or allow_all_taxa:
                raise ValueError("Taxon scoping applies only to UniProt idmapping")
            return write_knowledgebase(
                entries=self._required_entries(),
                canonical_sequences=self._canonical_sequences,
                isoform_sequences=self._isoform_sequences,
                release_version=self._release_version,
                path=Path(path),
                if_exists=if_exists,
            )
        self._require_mode("idmapping")
        normalized = normalize_taxids(tuple(taxon_ids or ()))
        if not normalized and not allow_all_taxa:
            raise ValueError("Writing all UniProt taxids requires allow_all_taxa=True")
        source = self._required_source_path()
        dataset = TidyDataset(
            resource_name="uniprot",
            frames={"mapping": self.scan_mapping(taxon_ids=normalized)},
            source=TidySource(
                "idmapping_selected",
                source,
                _mapping_media_type(self._mapping_kind or ""),
            ),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile=IDMAPPING_SOURCE_SCHEMA_PROFILE,
            build_id_prefix="uniprot-id-mapping",
            assets=(TidyAsset("mapping.parquet", "canonical", "mapping"),),
            release_version=self._release_version,
            extra_metadata={"bioextract.capability.mapping": "true"},
        )
        return dataset.write_duckdb(
            path,
            if_exists=if_exists,
            preserve_source_headers=("mapping",),
        )

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Open a caller-owned read-only DuckDB connection.

        Examples:
            >>> database.connect()  # doctest: +SKIP
            <duckdb.DuckDBPyConnection ...>
        """
        if self._mode not in {"knowledgebase_publication", "idmapping_publication"}:
            raise CapabilityError("connect() requires UniProtDatabase.from_duckdb()")
        if self._duckdb_path is None:
            raise RuntimeError("UniProt DuckDB path is unavailable")
        if _file_identity(self._duckdb_path) != self._duckdb_identity:
            raise IntegrityError(
                "UniProt publication was replaced; reopen it with from_duckdb()"
            )
        connection = duckdb.connect(str(self._duckdb_path), read_only=True)
        if _file_identity(self._duckdb_path) != self._duckdb_identity:
            connection.close()
            raise IntegrityError(
                "UniProt publication was replaced; reopen it with from_duckdb()"
            )
        return connection

    def select_ids(
        self,
        ids: Iterable[str],
        *,
        namespace: str,
        taxon_ids: Iterable[str | int] | None = None,
    ) -> UniProtSelection:
        """Select proteins through one supported identifier namespace.

        Examples:
            >>> database.select_ids(["P04637"], namespace="uniprot")  # doctest: +SKIP
            UniProtSelection(...)
        """
        self._require_mode("knowledgebase_publication")
        normalized_taxa = normalize_taxids(tuple(taxon_ids or ()))
        return make_selection(self, ids, namespace=namespace, taxon_ids=normalized_taxa)

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: str,
        taxon_ids: Iterable[str | int] | None = None,
    ) -> UniProtSelection:
        """Select identifiers while preserving caller group labels.

        Examples:
            >>> database.select_groups({"case": ["P04637"]}, namespace="uniprot")  # doctest: +SKIP
            UniProtSelection(...)
        """
        self._require_mode("knowledgebase_publication")
        normalized_taxa = normalize_taxids(tuple(taxon_ids or ()))
        return make_selection(
            self,
            (),
            namespace=namespace,
            groups=ids_by_group,
            taxon_ids=normalized_taxa,
        )

    def _scan_mapping(self, taxon_ids: tuple[str, ...]) -> pl.LazyFrame:
        if self._mode == "idmapping_publication":
            connection = self.connect()
            query = "SELECT * FROM mapping"
            params: list[str] = []
            if taxon_ids:
                query += f" WHERE TaxId IN ({','.join('?' for _ in taxon_ids)})"
                params.extend(taxon_ids)
            return connection.sql(query, params=params).pl(lazy=True)
        source = self._required_source_path()
        match self._mapping_kind:
            case "raw_plain" | "raw_gzip":
                return scan_raw_idmapping_selected(source)
            case "parquet":
                return pl.scan_parquet(source)
            case "hive":
                return scan_hive_mapping_dataset(source)
            case _:
                raise RuntimeError("Unknown UniProt idmapping kind")

    def _require_mode(self, mode: str) -> None:
        if self._mode != mode:
            raise CapabilityError(f"Operation requires UniProt {mode} capability")

    def _required_source_path(self) -> Path:
        if self._mode != "idmapping":
            raise CapabilityError("Operation requires UniProt idmapping source")
        if self._source_path is None:
            raise RuntimeError("UniProt idmapping source is unavailable")
        return self._source_path

    def _required_entries(self) -> Path:
        if self._entries is None:
            raise RuntimeError("UniProtKB entries source is unavailable")
        return self._entries


def _required_file(path: os.PathLike[str] | str, role: str) -> Path:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"UniProt {role} file not found: {source}")
    return source


def _normalize_version(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        raise ValueError("release_version must be non-empty when provided")
    return normalized


def _mapping_kind(path: Path) -> str:
    if path.is_dir():
        if not has_hive_parquet_candidates(path):
            raise ValueError(f"UniProt directory contains no parquet files: {path}")
        return "hive"
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic == b"PAR1":
        return "parquet"
    if magic.startswith(b"\x1f\x8b"):
        return "raw_gzip"
    return "raw_plain"


def _mapping_media_type(kind: str) -> str:
    return {
        "raw_plain": MEDIA_TYPE_TSV,
        "raw_gzip": MEDIA_TYPE_TSV_GZIP,
        "parquet": MEDIA_TYPE_PARQUET,
        "hive": MEDIA_TYPE_PARQUET_DATASET,
    }[kind]


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    try:
        stat = path.stat()
    except OSError as error:
        raise IntegrityError(
            "UniProt publication is unavailable; reopen it with from_duckdb()"
        ) from error
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns
