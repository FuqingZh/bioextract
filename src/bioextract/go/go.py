from __future__ import annotations

import json
import os
import tarfile
import zipfile
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal, cast

import duckdb
import polars as pl

from bioextract._publication import (
    BIOEXTRACT_RELATIONS,
    DuckDBWriteResult,
    validate_duckdb_metadata_v1,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError, IntegrityError

from .ontology.constant import (
    ASSET_SPECS,
    GO_SUBSET_GOSLIM_GENERIC,
    MEDIA_TYPE_OBO,
    SCHEMA_VERSION,
)
from .ontology.parse import (
    read_obo_subset_definitions,
    scan_obo_term_records,
    validate_go_id,
)
from .ontology.tidy import (
    _build_tidy_frames,  # pyright: ignore[reportPrivateUsage]  # owned ontology helper
    extract_subcell_frame,
)

__all__ = ["GODatabase"]


GoNamespace = Literal[
    "biological_process",
    "cellular_component",
    "molecular_function",
]


class _GoNamespace(StrEnum):
    BIOLOGICAL_PROCESS = "biological_process"
    CELLULAR_COMPONENT = "cellular_component"
    MOLECULAR_FUNCTION = "molecular_function"


class GoSubsetId(StrEnum):
    """Provide named GO subset IDs without restricting arbitrary subset text.

    `GODatabase.select_terms()` also accepts ordinary strings so snapshots can expose
    subsets added outside this convenience enum.

    Examples:
        Select the generic GO slim terms without spelling its raw subset ID:

        >>> db = GODatabase.from_obo("data/go-basic.obo")
        >>> db.select_terms(
        ...     subset_id=GoSubsetId.GOSLIM_GENERIC
        ... )["go_id"].to_list()
        ['GO:0000001', 'GO:0000002', 'GO:0005575', 'GO:0005737']
    """

    GOSLIM_GENERIC = GO_SUBSET_GOSLIM_GENERIC


@dataclass(frozen=True, slots=True)
class _GoSnapshot:
    file_obo: Path | None = None


@dataclass(slots=True)
class GODatabase:
    """Path-first access to a local Gene Ontology OBO snapshot.

    `GODatabase` is the public entrypoint for extracting tidy ontology tables from a
    local GO OBO file. It keeps the raw file path and builds materialized
    Polars frames only when an operation requests them.

    The default tidy output is a flat ontology snapshot with canonical term and
    edge tables plus derived graph tables. `extract_subcell()` is a convenience
    view over cellular component terms for subcellular-location workflows.

    Examples:
        Read a hierarchical edge from a compact local snapshot:

        >>> db = GODatabase.from_obo("data/go-basic.obo")
        >>> db.build_tidy().frames["edge"].filter(
        ...     pl.col("relation_type") == "is_a"
        ... ).select(
        ...     "child_go_id", "parent_go_id", "relation_type"
        ... ).head(1).collect().to_dicts()
        [{'child_go_id': 'GO:0000002', 'parent_go_id': 'GO:0000001', 'relation_type': 'is_a'}]
    """

    snapshot: _GoSnapshot
    _tidy: TidyDataset | None = field(default=None, init=False, repr=False)
    _publication_path: Path | None = field(default=None, init=False, repr=False)
    _publication_identity: tuple[int, int, int, int, int] | None = field(
        default=None, init=False, repr=False
    )

    @classmethod
    def from_obo(
        cls,
        path: os.PathLike[str] | str,
    ) -> GODatabase:
        """Create a dataset handle from a local GO OBO file.

        Args:
            path: Local Gene Ontology OBO path or supported archive.

        Returns:
            A dataset handle that can build tidy ontology frames and subcellular
            component exports.

        Raises:
            FileNotFoundError: If the OBO file does not exist.

        Examples:
            Open a local GO fixture and read one parsed term:

            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> db.select_terms(term_ids=["GO:0000002"])["term_name"].item()
            'child process'
        """
        source_path = _resolve_obo_input(Path(path))
        return cls(snapshot=_GoSnapshot(file_obo=source_path))

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> GODatabase:
        """Open a validated GO ontology publication for domain and SQL access.

        Args:
            path: A bioextract GO metadata-v1 DuckDB publication.

        Returns:
            A publication-backed handle pinned to the validated file identity.

        Raises:
            FileNotFoundError: If the publication does not exist.
            IntegrityError: If the publication contract is invalid or the file
                changes while it is being validated.

        Examples:
            Reopen a publication for ontology selection:

            >>> db = GODatabase.from_duckdb("tidy/go.duckdb")  # doctest: +SKIP
            >>> db.select_terms(term_ids=["GO:0005575"]).height  # doctest: +SKIP
            1
        """
        publication_path = Path(path).resolve()
        identity_before = _file_identity(publication_path)
        try:
            _validate_go_publication(publication_path)
        except (KeyError, TypeError, ValueError) as error:
            raise IntegrityError(str(error)) from error
        identity_after = _file_identity(publication_path)
        if identity_after != identity_before:
            raise IntegrityError("GO publication changed during validation")
        result = cls(snapshot=_GoSnapshot())
        result._publication_path = publication_path
        result._publication_identity = identity_after
        return result

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a fresh caller-owned read-only DuckDB connection.

        Raises:
            CapabilityError: If this handle was created from an OBO source.
            IntegrityError: If the validated publication path was replaced.

        Examples:
            Query the publication through native DuckDB SQL:

            >>> db = GODatabase.from_duckdb("tidy/go.duckdb")  # doctest: +SKIP
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM term").fetchone()[0]
            >>> count >= 0  # doctest: +SKIP
            True
        """
        path = self._publication_path
        if path is None:
            raise CapabilityError("connect() requires GODatabase.from_duckdb()")
        self._assert_publication_identity()
        connection = duckdb.connect(str(path), read_only=True)
        try:
            self._assert_publication_identity()
        except BaseException:
            connection.close()
            raise
        return connection

    def build_tidy(self) -> TidyDataset:
        """Build the GO tidy dataset.

        Returns:
            A `TidyDataset` with `term`, `edge`, `synonym`, `xref`, `alt_id`,
            `subset_membership`, `subset_definition`, `ancestor_all`, and
            `depth` frames.

        Examples:
            Inspect the frame names built from a local OBO snapshot:

            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> sorted(db.build_tidy().frames)
            ['alt_id', 'ancestor_all', 'depth', 'edge', 'subset_definition', 'subset_membership', 'synonym', 'term', 'xref']
        """
        if self._tidy is not None:
            return self._tidy

        if self._publication_path is not None:
            frames = {
                frame_name: frame.lazy()
                for frame_name, frame in self._read_publication_frames().items()
            }
            self._tidy = TidyDataset(
                frames=frames,
                source=(),
                resource_schema_version=SCHEMA_VERSION,
                source_schema_profile="gene-ontology-obo-v1",
                build_id_prefix="go-ontology-publication",
                assets=tuple(
                    TidyAsset(path=path, kind=kind, frame_name=frame_name)
                    for path, kind, frame_name in ASSET_SPECS
                ),
                resource_name="go",
            )
            return self._tidy

        if self.snapshot.file_obo is None:
            raise CapabilityError("GO OBO source is unavailable")
        records = scan_obo_term_records(self.snapshot.file_obo)
        subset_definitions = read_obo_subset_definitions(self.snapshot.file_obo)
        frames = {
            frame_name: frame.lazy()
            for frame_name, frame in _build_tidy_frames(
                records,
                subset_definitions=subset_definitions,
            ).items()
        }
        self._tidy = TidyDataset(
            frames=frames,
            source=TidySource(
                logical_name="go_obo",
                path=self.snapshot.file_obo,
                media_type=_obo_media_type(self.snapshot.file_obo),
            ),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="gene-ontology-obo-v1",
            build_id_prefix=f"go-ontology-{self.snapshot.file_obo.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
            resource_name="go",
        )
        return self._tidy

    def select_terms(
        self,
        *,
        term_ids: Iterable[str] | None = None,
        namespace: GoNamespace | None = None,
        subset_id: str | GoSubsetId | None = None,
        include_obsolete: bool = False,
        resolve_alt_ids: bool = True,
    ) -> pl.DataFrame:
        """Select GO terms from the current OBO snapshot.

        Args:
            term_ids: Optional GO IDs to select. Alternate GO IDs are resolved
                to primary GO IDs by default.
            namespace: Optional GO namespace filter.
            subset_id: Optional OBO subset membership filter.
            include_obsolete: Whether to keep obsolete GO terms.
            resolve_alt_ids: Whether to resolve alternate GO IDs through
                the `alt_id` frame.

        Returns:
            A stable term view. When `term_ids` is provided, `input_go_id` is
            included to show canonicalization.

        Examples:
            Resolve an alternate ID from the compact GO fixture:

            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> db.select_terms(term_ids=["GO:1234567"]).select(
            ...     "input_go_id", "go_id"
            ... ).to_dicts()
            [{'input_go_id': 'GO:1234567', 'go_id': 'GO:0000002'}]
        """
        frame_names = {"term", "alt_id"}
        if subset_id is not None:
            frame_names.add("subset_membership")
        frames = self._collect_frames(frame_names)
        df_term = frames["term"]
        if not include_obsolete:
            df_term = df_term.filter(~pl.col("is_obsolete"))

        is_term_id_filter = term_ids is not None
        if term_ids is not None:
            df_input_terms = create_go_term_input_frame(term_ids)
            df_term = select_terms_by_ids(
                df_term,
                frames["alt_id"],
                df_input_terms,
                resolve_alt_ids=resolve_alt_ids,
            )

        if namespace is not None:
            namespace_value = normalize_go_namespace(namespace)
            df_term = df_term.filter(pl.col("namespace") == namespace_value)

        subset_id_value = None if subset_id is None else normalize_subset_id(subset_id)
        if subset_id_value is not None:
            df_subset = frames["subset_membership"].filter(
                pl.col("subset_id") == subset_id_value
            )
            df_term = df_term.join(df_subset, on="go_id", how="inner")

        cols_out = [
            "go_id",
            "term_name",
            "namespace",
            "definition",
            "is_obsolete",
            "comment",
        ]
        if is_term_id_filter:
            cols_out.insert(0, "input_go_id")
        if subset_id_value is not None:
            cols_out.append("subset_id")

        df_selected = df_term.select(cols_out)
        if is_term_id_filter and "input_order" in df_term.columns:
            return (
                df_term.select(cols_out + ["input_order"])
                .sort("input_order", "go_id")
                .drop("input_order")
            )
        return df_selected.sort("go_id")

    def list_subsets(self) -> pl.DataFrame:
        """List OBO subset definitions and term counts in this snapshot.

        Returns:
            A table sorted by `subset_id` with `subset_id`, `subset_name`, and
            `num_terms`. Declared subsets with no members are retained with a
            zero count.

        Examples:
            Discover a subset's display name and term count:

            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> db.list_subsets().row(0, named=True)
            {'subset_id': 'goslim_generic', 'subset_name': 'Generic GO slim', 'num_terms': 5}
        """
        frames = self._collect_frames({"subset_membership", "subset_definition"})
        df_membership = frames["subset_membership"]
        df_definition = frames["subset_definition"]
        df_counts = (
            df_membership.group_by("subset_id").agg(
                pl.col("go_id").n_unique().alias("num_terms")
            )
            if df_membership.height
            else pl.DataFrame(
                {
                    "subset_id": pl.Series([], dtype=pl.String),
                    "num_terms": pl.Series([], dtype=pl.UInt32),
                }
            )
        )
        df_subsets = df_definition.join(df_counts, on="subset_id", how="full").select(
            pl.coalesce("subset_id", "subset_id_right").alias("subset_id"),
            "subset_name",
            pl.col("num_terms").fill_null(0),
        )
        return df_subsets.sort("subset_id")

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: Literal["fail", "replace"] = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish the complete ontology as one DuckDB database.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(Path(dir_out) / "go.duckdb")
            ...     "term_relation" in result.tables
            True
        """
        if self._publication_path is not None:
            raise CapabilityError("write_duckdb() requires a GO OBO source handle")
        return self.build_tidy().write_duckdb(
            Path(path),
            table_names={
                "edge": "term_relation",
                "synonym": "term_synonym",
                "xref": "term_xref",
                "alt_id": "term_alternate_id",
                "ancestor_all": "term_ancestor",
                "depth": "term_depth",
            },
            if_exists=if_exists,
        )

    def extract_subcell(self, *, include_obsolete: bool = False) -> pl.DataFrame:
        """Extract non-obsolete cellular component terms as a subcell table.

        Args:
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            A DataFrame with GO ID, subcell name, definition, and depth columns.

        Examples:
            Extract cellular-component IDs from the compact fixture:

            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> db.extract_subcell()["go_id"].to_list()
            ['GO:0005575', 'GO:0005737']
        """
        return extract_subcell_frame(
            self._collect_frames({"term", "depth"}),
            include_obsolete=include_obsolete,
        )

    def write_subcell(
        self,
        path: os.PathLike[str] | str,
        *,
        include_obsolete: bool = False,
    ) -> Path:
        """Write the cellular component subcell table as a parquet file.

        Args:
            path: Output parquet path.
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            The output path that was written.

        Examples:
            Write the subcell projection and read back its public columns:

            >>> from tempfile import TemporaryDirectory
            >>> db = GODatabase.from_obo("data/go-basic.obo")
            >>> with TemporaryDirectory() as dir_out:
            ...     path = Path(dir_out) / "subcell.parquet"
            ...     _ = db.write_subcell(path)
            ...     pl.read_parquet(path).select(
            ...         "go_id", "subcell_name"
            ...     ).to_dicts()
            [{'go_id': 'GO:0005575', 'subcell_name': 'cellular_component'}, {'go_id': 'GO:0005737', 'subcell_name': 'cytoplasm'}]
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.extract_subcell(include_obsolete=include_obsolete).lazy().sink_parquet(
            path
        )
        return path

    def _assert_publication_identity(self) -> None:
        path = self._publication_path
        if path is None or _file_identity(path) != self._publication_identity:
            raise IntegrityError(
                "GO publication was replaced; reopen it with from_duckdb()"
            )

    def _collect_frames(self, frame_names: set[str]) -> dict[str, pl.DataFrame]:
        if self._publication_path is not None:
            return self._read_publication_frames(frame_names)
        return {
            frame_name: self.build_tidy().frames[frame_name].collect()
            for frame_name in frame_names
        }

    def _read_publication_frames(
        self, frame_names: set[str] | None = None
    ) -> dict[str, pl.DataFrame]:
        table_names = {
            "term": "term",
            "edge": "term_relation",
            "synonym": "term_synonym",
            "xref": "term_xref",
            "alt_id": "term_alternate_id",
            "subset_membership": "subset_membership",
            "subset_definition": "subset_definition",
            "ancestor_all": "term_ancestor",
            "depth": "term_depth",
        }
        selected_names = set(table_names) if frame_names is None else frame_names
        with self.connect() as connection:
            return {
                frame_name: pl.read_database(  # pyright: ignore[reportUnknownMemberType]
                    f'SELECT * FROM "{table_names[frame_name]}"', connection
                )
                for frame_name in selected_names
            }


def normalize_go_namespace(namespace: str) -> str:
    try:
        return _GoNamespace(namespace).value
    except ValueError as error:
        values = ", ".join(namespace.value for namespace in _GoNamespace)
        raise ValueError(
            f"Unsupported GO namespace: {namespace!r}; expected {values}"
        ) from error


def normalize_subset_id(subset_id: str | GoSubsetId) -> str:
    return str(
        subset_id.value if isinstance(subset_id, GoSubsetId) else subset_id
    ).strip()


def _resolve_obo_input(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"GO OBO file not found: {path}")
    if path.is_file():
        return path
    candidates = [
        candidate
        for candidate in path.rglob("*")
        if candidate.is_file()
        and (
            ".obo" in candidate.name.lower()
            or candidate.name.lower().endswith((".zip", ".tar", ".tgz"))
        )
    ]
    if len(candidates) != 1:
        raise ValueError(
            "GO ontology directory must contain exactly one recognizable "
            f"ontology input; found {len(candidates)}"
        )
    return candidates[0]


def _obo_media_type(path: Path) -> str:
    with path.open("rb") as handle:
        if handle.read(2) == b"\x1f\x8b":
            return "application/gzip"
    if zipfile.is_zipfile(path):
        return "application/zip"
    if tarfile.is_tarfile(path):
        return "application/x-tar"
    return MEDIA_TYPE_OBO


def create_go_term_input_frame(term_ids: Iterable[str]) -> pl.DataFrame:
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for term_id in term_ids:
        input_go_id = str(term_id).strip()
        if not input_go_id or input_go_id in seen:
            continue
        validate_go_id(input_go_id)
        seen.add(input_go_id)
        rows.append({"input_go_id": input_go_id, "input_order": len(rows)})
    return pl.DataFrame(
        rows,
        schema={
            "input_go_id": pl.String,
            "input_order": pl.UInt32,
        },
    )


def select_terms_by_ids(
    df_term: pl.DataFrame,
    df_alt_id: pl.DataFrame,
    df_input_terms: pl.DataFrame,
    *,
    resolve_alt_ids: bool,
) -> pl.DataFrame:
    if df_input_terms.height == 0:
        return df_term.head(0).with_columns(
            pl.lit(None, dtype=pl.String).alias("input_go_id"),
            pl.lit(None, dtype=pl.UInt32).alias("input_order"),
        )

    df_primary = df_input_terms.rename({"input_go_id": "go_id"}).join(
        df_term.select("go_id"),
        on="go_id",
        how="inner",
    )
    df_primary = df_primary.with_columns(pl.col("go_id").alias("input_go_id"))

    if resolve_alt_ids:
        df_alt = (
            df_input_terms.join(
                df_alt_id,
                left_on="input_go_id",
                right_on="alt_go_id",
                how="inner",
            )
            .rename({"primary_go_id": "go_id"})
            .select("input_go_id", "input_order", "go_id")
        )
    else:
        df_alt = pl.DataFrame(
            {
                "input_go_id": pl.Series([], dtype=pl.String),
                "input_order": pl.Series([], dtype=pl.UInt32),
                "go_id": pl.Series([], dtype=pl.String),
            }
        )

    df_term_ids = pl.concat(
        [
            df_primary.select("input_go_id", "input_order", "go_id"),
            df_alt,
        ],
        how="vertical",
    ).unique(subset=["input_go_id", "go_id"], keep="first", maintain_order=True)
    return df_term_ids.join(df_term, on="go_id", how="inner")


_GO_TABLE_CONTRACTS: dict[str, tuple[str, tuple[tuple[str, str], ...]]] = {
    "term": (
        "canonical",
        (
            ("go_id", "VARCHAR"),
            ("term_name", "VARCHAR"),
            ("namespace", "VARCHAR"),
            ("definition", "VARCHAR"),
            ("is_obsolete", "BOOLEAN"),
            ("comment", "VARCHAR"),
        ),
    ),
    "term_relation": (
        "canonical",
        (
            ("child_go_id", "VARCHAR"),
            ("parent_go_id", "VARCHAR"),
            ("relation_type", "VARCHAR"),
            ("source_clause", "VARCHAR"),
        ),
    ),
    "term_synonym": (
        "canonical",
        (
            ("go_id", "VARCHAR"),
            ("synonym_text", "VARCHAR"),
            ("synonym_scope", "VARCHAR"),
            ("synonym_type_name", "VARCHAR"),
            ("dbxref_text", "VARCHAR"),
        ),
    ),
    "term_xref": (
        "canonical",
        (
            ("go_id", "VARCHAR"),
            ("xref_text", "VARCHAR"),
            ("xref_db", "VARCHAR"),
            ("xref_id", "VARCHAR"),
        ),
    ),
    "term_alternate_id": (
        "canonical",
        (("alt_go_id", "VARCHAR"), ("primary_go_id", "VARCHAR")),
    ),
    "subset_membership": (
        "canonical",
        (("go_id", "VARCHAR"), ("subset_id", "VARCHAR")),
    ),
    "subset_definition": (
        "canonical",
        (("subset_id", "VARCHAR"), ("subset_name", "VARCHAR")),
    ),
    "term_ancestor": (
        "derived",
        (
            ("go_id", "VARCHAR"),
            ("ancestor_go_id", "VARCHAR"),
            ("min_distance", "INTEGER"),
        ),
    ),
    "term_depth": (
        "derived",
        (
            ("go_id", "VARCHAR"),
            ("namespace", "VARCHAR"),
            ("min_depth_from_root", "INTEGER"),
            ("max_depth_from_root", "INTEGER"),
        ),
    ),
}


def _validate_go_publication(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with duckdb.connect(str(path), read_only=True) as connection:
            metadata_rows = connection.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
            metadata = {str(row[0]): str(row[1]) for row in metadata_rows}
            if len(metadata) != len(metadata_rows):
                raise ValueError("GO publication has duplicate metadata keys")
            if metadata.get("bioextract.metadata_schema_version") != "1":
                raise ValueError("Unsupported GO metadata schema version")
            validate_duckdb_metadata_v1(connection, metadata)
            if metadata.get("bioextract.resource_name") != "go":
                raise ValueError("DuckDB file is not a bioextract GO publication")
            if metadata.get("bioextract.source_schema_profile") != (
                "gene-ontology-obo-v1"
            ):
                raise ValueError("Unsupported GO source schema profile")
            if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported GO resource schema version")

            source_rows = connection.execute(
                "SELECT logical_name, display_path, bytes, media_type, sha256 "
                "FROM _bioextract.source_file"
            ).fetchall()
            if len(source_rows) != 1 or source_rows[0][0] != "go_obo":
                raise ValueError("GO source role inventory is unsupported")
            embedded_sources: object = json.loads(metadata["bioextract.sources"])
            if (
                not isinstance(embedded_sources, list)
                or len(cast(list[object], embedded_sources)) != 1
            ):
                raise ValueError("GO embedded source inventory is unsupported")
            if int(source_rows[0][2]) < 0:
                raise ValueError("GO source byte count is unsupported")

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
            } | {("main", name, "BASE TABLE") for name in _GO_TABLE_CONTRACTS}
            if relations != expected_relations:
                raise ValueError("GO physical table/view inventory is unsupported")

            info_rows = connection.execute(
                "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
            ).fetchall()
            recorded = {str(row[0]): (str(row[1]), int(row[2])) for row in info_rows}
            if len(recorded) != len(info_rows) or set(recorded) != set(
                _GO_TABLE_CONTRACTS
            ):
                raise ValueError("GO table inventory does not match metadata")
            for table_name, (role, row_count) in recorded.items():
                expected_role, expected_schema = _GO_TABLE_CONTRACTS[table_name]
                if role != expected_role or row_count < 0:
                    raise ValueError("GO table capability inventory is unsupported")
                actual_schema = tuple(
                    (str(row[1]), str(row[2]))
                    for row in connection.execute(
                        f'PRAGMA table_info("{table_name}")'
                    ).fetchall()
                )
                if actual_schema != expected_schema:
                    raise ValueError(f"GO table schema is unsupported: {table_name}")
            if connection.execute(
                "SELECT count(*) FROM _bioextract.column_mapping"
            ).fetchone() != (0,):
                raise ValueError("GO column provenance inventory is unsupported")
    except duckdb.Error as error:
        raise ValueError(f"Cannot open GO DuckDB publication: {path}") from error


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
