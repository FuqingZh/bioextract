from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal

import polars as pl

from bioextract._shared import validate_file_size
from bioextract._tidy import TidyAsset, TidyDataset, TidySource, TidyWriteReport

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
from .ontology.tidy import _build_tidy_frames, extract_subcell_frame

__all__ = [
    "GoDb",
    "GoNamespace",
    "GoResourceLimits",
    "GoSubsetId",
    "GoTidyDataset",
]


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

    `GoDb.select_terms()` also accepts ordinary strings so snapshots can expose
    subsets added outside this convenience enum.

    Examples:
        Select the generic GO slim terms without spelling its raw subset ID:

        >>> db = GoDb.from_obo("data/go-basic.obo")
        >>> db.select_terms(
        ...     subset_id=GoSubsetId.GOSLIM_GENERIC
        ... )["go_id"].to_list()
        ['GO:0000001', 'GO:0000002', 'GO:0005575', 'GO:0005737']
    """

    GOSLIM_GENERIC = GO_SUBSET_GOSLIM_GENERIC


@dataclass(frozen=True, slots=True)
class GoResourceLimits:
    """Configure fail-fast limits for a GO OBO snapshot.

    Attributes:
        file_obo_bytes_max: Maximum on-disk OBO file size in bytes, or `None`
            to disable the size check.

    Examples:
        Reject a snapshot before parsing when it exceeds the configured limit:

        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_obo = Path(dir_tmp) / "go-basic.obo"
        ...     _ = file_obo.write_text("format-version: 1.2\\n", encoding="utf-8")
        ...     limits = GoResourceLimits(file_obo_bytes_max=1)
        ...     try:
        ...         GoDb.from_obo(file_obo, limits=limits)
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_obo_bytes_max: int | None = None


@dataclass(frozen=True, slots=True)
class _GoSnapshot:
    file_obo: Path


GoTidyDataset = TidyDataset


@dataclass(slots=True)
class GoDb:
    """Path-first access to a local Gene Ontology OBO snapshot.

    `GoDb` is the public entrypoint for extracting tidy ontology tables from a
    local GO OBO file. It keeps the raw file path and resource limits, then
    builds materialized Polars frames only when tidy or convenience exports are
    requested.

    The default tidy output is a flat ontology snapshot with canonical term and
    edge tables plus derived graph tables. `extract_subcell()` is a convenience
    view over cellular component terms for subcellular-location workflows.

    Examples:
        Read a hierarchical edge from a compact local snapshot:

        >>> db = GoDb.from_obo("data/go-basic.obo")
        >>> db.build_tidy().frames["edge"].filter(
        ...     pl.col("relation_type") == "is_a"
        ... ).select(
        ...     "child_go_id", "parent_go_id", "relation_type"
        ... ).head(1).collect().to_dicts()
        [{'child_go_id': 'GO:0000002', 'parent_go_id': 'GO:0000001', 'relation_type': 'is_a'}]
    """

    snapshot: _GoSnapshot
    limits: GoResourceLimits
    _tidy: GoTidyDataset | None = field(default=None, init=False, repr=False)

    DEFAULT_RESOURCE_LIMITS = GoResourceLimits()

    @classmethod
    def from_obo(
        cls,
        file_obo: os.PathLike[str] | str,
        *,
        limits: GoResourceLimits | None = None,
    ) -> GoDb:
        """Create a dataset handle from a local GO OBO file.

        Args:
            file_obo: Path to a local Gene Ontology OBO file.
            limits: Optional resource policy. When omitted, the file-size
                check is disabled.

        Returns:
            A dataset handle that can build tidy ontology frames and subcellular
            component exports.

        Raises:
            FileNotFoundError: If the OBO file does not exist.
            ValueError: If the configured file-size limit is exceeded.

        Examples:
            Open a local GO fixture and read one parsed term:

            >>> db = GoDb.from_obo("data/go-basic.obo")
            >>> db.select_terms(term_ids=["GO:0000002"])["term_name"].item()
            'child process'
        """
        file_obo = Path(file_obo)
        if not file_obo.exists():
            raise FileNotFoundError(f"GO OBO file not found: {file_obo}")

        limits_resolved = GoResourceLimits() if limits is None else limits
        validate_file_size(
            file_path=file_obo,
            size_max=limits_resolved.file_obo_bytes_max,
            label="GO OBO file",
        )
        return cls(snapshot=_GoSnapshot(file_obo=file_obo), limits=limits_resolved)

    def build_tidy(self) -> GoTidyDataset:
        """Build the GO tidy dataset.

        Returns:
            A `TidyDataset` with `term`, `edge`, `synonym`, `xref`, `alt_id`,
            `subset_membership`, `subset_definition`, `ancestor_all`, and
            `depth` frames.

        Examples:
            Inspect the frame names built from a local OBO snapshot:

            >>> db = GoDb.from_obo("data/go-basic.obo")
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
        self._tidy = GoTidyDataset(
            frames=frames,
            source=TidySource(path=self.snapshot.file_obo, media_type=MEDIA_TYPE_OBO),
            schema_version=SCHEMA_VERSION,
            build_id_prefix=f"go-ontology-{self.snapshot.file_obo.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
        )
        return self._tidy

    def select_terms(
        self,
        *,
        term_ids: Iterable[str] | None = None,
        namespace: GoNamespace | None = None,
        subset_id: str | GoSubsetId | None = None,
        include_obsolete: bool = False,
        should_resolve_alt_ids: bool = True,
    ) -> pl.DataFrame:
        """Select GO terms from the current OBO snapshot.

        Args:
            term_ids: Optional GO IDs to select. Alternate GO IDs are resolved
                to primary GO IDs by default.
            namespace: Optional GO namespace filter.
            subset_id: Optional OBO subset membership filter.
            include_obsolete: Whether to keep obsolete GO terms.
            should_resolve_alt_ids: Whether to resolve alternate GO IDs through
                the `alt_id` frame.

        Returns:
            A stable term view. When `term_ids` is provided, `input_go_id` is
            included to show canonicalization.

        Examples:
            Resolve an alternate ID from the compact GO fixture:

            >>> db = GoDb.from_obo("data/go-basic.obo")
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
                should_resolve_alt_ids=should_resolve_alt_ids,
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

            >>> db = GoDb.from_obo("data/go-basic.obo")
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

    def write_tidy(
        self,
        dir_out: os.PathLike[str] | str,
        *,
        should_write_manifest: bool = False,
        should_hash_assets: bool = False,
    ) -> TidyWriteReport:
        """Write the GO tidy dataset as flat parquet files.

        Args:
            dir_out: Output directory for parquet assets.
            should_write_manifest: Whether to write `manifest.json`.
            should_hash_assets: Whether to calculate asset checksums in the
                manifest.

        Returns:
            A write report with asset paths and optional manifest content.

        Examples:
            Write the nine declared GO assets:

            >>> db = GoDb.from_obo("data/go-basic.obo")
            >>> report = db.write_tidy("build/go-basic")
            >>> (report.assets[0].path, report.assets[-1].path)
            ('term.parquet', 'depth.parquet')
        """
        return self.build_tidy().write(
            Path(dir_out),
            should_write_manifest=should_write_manifest,
            should_hash_assets=should_hash_assets,
        )

    def extract_subcell(self, *, include_obsolete: bool = False) -> pl.DataFrame:
        """Extract non-obsolete cellular component terms as a subcell table.

        Args:
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            A DataFrame with GO ID, subcell name, definition, and depth columns.

        Examples:
            Extract cellular-component IDs from the compact fixture:

            >>> db = GoDb.from_obo("data/go-basic.obo")
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
        file_out: os.PathLike[str] | str,
        *,
        include_obsolete: bool = False,
    ) -> Path:
        """Write the cellular component subcell table as a parquet file.

        Args:
            file_out: Output parquet path.
            include_obsolete: Whether to keep obsolete cellular component terms.

        Returns:
            The output path that was written.

        Examples:
            Write the subcell projection and read back its public columns:

            >>> from tempfile import TemporaryDirectory
            >>> db = GoDb.from_obo("data/go-basic.obo")
            >>> with TemporaryDirectory() as dir_out:
            ...     file_out = Path(dir_out) / "subcell.parquet"
            ...     _ = db.write_subcell(file_out)
            ...     pl.read_parquet(file_out).select(
            ...         "go_id", "subcell_name"
            ...     ).to_dicts()
            [{'go_id': 'GO:0005575', 'subcell_name': 'cellular_component'}, {'go_id': 'GO:0005737', 'subcell_name': 'cytoplasm'}]
        """
        file_out = Path(file_out)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        self.extract_subcell(include_obsolete=include_obsolete).lazy().sink_parquet(
            file_out
        )
        return file_out


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
    should_resolve_alt_ids: bool,
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

    if should_resolve_alt_ids:
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
