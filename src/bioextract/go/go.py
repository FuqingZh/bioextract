from __future__ import annotations

import os
import tarfile
import zipfile
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal

import polars as pl

from bioextract._publication import DuckDBWriteResult
from bioextract._tidy import TidyAsset, TidyDataset, TidySource

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
    file_obo: Path


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
        frames = {
            frame_name: frame.collect()
            for frame_name, frame in self.build_tidy().frames.items()
        }
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
        frames = {
            frame_name: frame.collect()
            for frame_name, frame in self.build_tidy().frames.items()
        }
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
            {
                frame_name: frame.collect()
                for frame_name, frame in self.build_tidy().frames.items()
            },
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
