from __future__ import annotations

import copy
import json
import os
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import (
    _RelationScanRequest,  # pyright: ignore[reportPrivateUsage]  # typed source boundary
    register_replayable_source,
)
from bioextract._publication import (
    BIOEXTRACT_RELATIONS,
    METADATA_SCHEMA_VERSION,
    DuckDBWriteResult,
    validate_duckdb_metadata_v2,
)
from bioextract._shared import (
    create_group_input_frames,
    create_input_id_frame,
)
from bioextract._tidy import TidyAsset, TidyDataset, TidySource
from bioextract.errors import CapabilityError, IntegrityError

from .constant import (
    ASSET_SPECS,
    MEDIA_TYPE_GMT,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUP_MAPPING,
    SCHEMA_GROUPS,
    SCHEMA_MAPPING,
    SCHEMA_PATHWAY,
    SCHEMA_TERM2GENE,
    SCHEMA_TERM2NAME,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    iter_gmt_records,
    read_gmt_frames,
    resolve_gmt_sources,
)

__all__ = [
    "WikiPathwaysDatabase",
]


@dataclass(frozen=True, slots=True)
class _WikiPathwaysSnapshot:
    files_gmt: tuple[Path, ...]
    species: str | None = None


class _ReopenedWikiPathwaysTidyDataset(TidyDataset):
    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        table_names: Mapping[str, str] | None = None,
        if_exists: str = "fail",
        source_columns: Mapping[str, Collection[str]] | None = None,
    ) -> DuckDBWriteResult:
        del path, table_names, if_exists, source_columns
        raise CapabilityError(
            "write_duckdb() requires a WikiPathways GMT source handle"
        )


@dataclass(slots=True)
class WikiPathwaysDatabase:
    """Source-first access to a local WikiPathways GMT snapshot.

    `WikiPathwaysDatabase` is the public entrypoint for extracting pathway gene sets
    from one or more local WikiPathways GMT files. GMT rows are interpreted as
    pathway metadata followed by NCBI Entrez Gene IDs; the class does not
    perform identifier conversion or enrichment statistics.

    Construct instances with :meth:`from_gmt`, then either extract whole
    resource frames or create single/grouped Entrez ID selections through
    :meth:`select_ids` and :meth:`select_groups`.

    Examples:
        Extract pathway mappings while preserving unmapped Entrez IDs:

        >>> db = WikiPathwaysDatabase.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
        ... ).with_species("Homo sapiens")
        >>> selection = db.select_ids(["2687", "MISSING"])
        >>> (
        ...     selection.mappings().collect()
        ...     .select("input_id", "wiki_pathways_id")
        ...     .to_dicts()
        ... )
        [{'input_id': '2687', 'wiki_pathways_id': 'WP100'}]
        >>> selection.unmatched_ids().collect().to_dicts()
        [{'input_id': 'MISSING'}]
    """

    snapshot: _WikiPathwaysSnapshot
    _frames: dict[str, pl.LazyFrame] | None = field(
        default=None, init=False, repr=False
    )
    _release_version: str | None = field(default=None, init=False, repr=False)
    _publication_path: Path | None = field(default=None, init=False, repr=False)
    _publication_identity: tuple[int, int, int, int, int] | None = field(
        default=None, init=False, repr=False
    )

    @classmethod
    def from_gmt(
        cls,
        source: os.PathLike[str] | str | Sequence[os.PathLike[str] | str],
        *,
        glob: bool = True,
    ) -> WikiPathwaysDatabase:
        """Create a dataset handle from one or more local WikiPathways GMT files.

        Args:
            source: One local path, a sequence of local paths, or glob
                expressions. Every resolved physical file must be unique.
            glob: Expand each source entry as a glob expression. If false,
                treat every entry literally.

        Returns:
            A dataset handle that can build pathway, term2gene, and term2name
            frames.

        Raises:
            FileNotFoundError: If a literal file does not exist or a glob
                expression has no matches.
            ValueError: At construction, if the source is empty, resolves to a
                directory, or repeats a physical file. GMT content parsing is deferred;
                malformed content, inconsistent Collection or Version values,
                and duplicate WikiPathways IDs raise when a frame is first
                extracted, built, or written.

        Examples:
            Open all GMT files in one local snapshot and retain human rows:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/*.gmt",
            ... ).with_species("Homo sapiens")
            >>> (
            ...     db.pathways().collect()
            ...     .select("wiki_pathways_id", "pathway_name")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'wiki_pathways_id': 'WP100', 'pathway_name': 'Glutathione metabolism'}]
        """
        files_gmt = resolve_gmt_sources(source, glob=glob)
        return cls(snapshot=_WikiPathwaysSnapshot(files_gmt=files_gmt))

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> WikiPathwaysDatabase:
        """Open a validated WikiPathways publication for domain and SQL access.

        Validation reads bounded metadata and catalog schemas only. The handle
        is pinned to the exact file that passed validation.

        Args:
            path: A bioextract WikiPathways metadata-v2 DuckDB publication.

        Returns:
            A publication-backed handle for extraction and selection.

        Raises:
            FileNotFoundError: If the publication does not exist.
            IntegrityError: If its metadata, inventory, or physical schema is
                incompatible, or if the file changes during validation.

        Examples:
            Reopen a publication and extract pathway metadata:

            >>> db = WikiPathwaysDatabase.from_duckdb(  # doctest: +SKIP
            ...     "tidy/wikipathways.duckdb"
            ... )
            >>> db.pathways().collect().height > 0  # doctest: +SKIP
            True
        """
        publication_path = Path(path).absolute()
        identity_before = _file_identity(publication_path)
        try:
            release_version, declared_species = _validate_wikipathways_publication(
                publication_path
            )
            identity_after = _file_identity(publication_path)
        except (KeyError, TypeError, ValueError) as error:
            raise IntegrityError(str(error)) from error
        except OSError as error:
            raise IntegrityError(
                "WikiPathways publication changed during validation"
            ) from error
        if identity_after != identity_before:
            raise IntegrityError("WikiPathways publication changed during validation")
        result = cls(
            snapshot=_WikiPathwaysSnapshot(
                files_gmt=(),
                species=declared_species,
            )
        )
        result._release_version = release_version
        result._publication_path = publication_path
        result._publication_identity = identity_after
        return result

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a fresh caller-owned read-only DuckDB connection.

        Raises:
            CapabilityError: If this handle was created from GMT source files.
            IntegrityError: If the validated publication was replaced or became
                unavailable.

        Examples:
            Run native SQL against a reopened publication:

            >>> db = WikiPathwaysDatabase.from_duckdb(  # doctest: +SKIP
            ...     "tidy/wikipathways.duckdb"
            ... )
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
            >>> count >= 0  # doctest: +SKIP
            True
        """
        path = self._publication_path
        if path is None:
            raise CapabilityError(
                "connect() requires WikiPathwaysDatabase.from_duckdb()"
            )
        self._assert_publication_identity()
        try:
            connection = duckdb.connect(str(path), read_only=True)
        except duckdb.Error as error:
            raise IntegrityError(
                "WikiPathways publication became unavailable; reopen it with "
                "from_duckdb()"
            ) from error
        try:
            self._assert_publication_identity()
        except BaseException:
            connection.close()
            raise
        return connection

    def select_ids(self, ids: Iterable[str]) -> WikiPathwaysSelection:
        """Create a single-query selection from Entrez Gene IDs.

        Args:
            ids: Input Entrez Gene IDs. Values are normalized with the shared
                input-ID normalizer.

        Returns:
            A selection that can extract pathway mappings and unmapped IDs.

        Examples:
            Select the pathway containing Entrez Gene ID 2687:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.select_ids(["2687"])
            ...     .mappings().collect()
            ...     .select("input_id", "wiki_pathways_id")
            ...     .to_dicts()
            ... )
            [{'input_id': '2687', 'wiki_pathways_id': 'WP100'}]
        """
        self._assert_publication_current()
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return WikiPathwaysSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
        )

    @property
    def species(self) -> str | None:
        """Return the exact species scope of this resource view.

        Examples:
            >>> db.species  # doctest: +SKIP
            'Homo sapiens'
        """

        return self.snapshot.species

    def with_species(self, species: str) -> WikiPathwaysDatabase:
        """Create a species-scoped view sharing this snapshot identity.

        The scope is applied to ``pathway`` first and then semi-joins
        ``pathway_gene`` by pathway ID. This prevents an Entrez Gene ID shared
        by two species from carrying a pathway across species boundaries.

        Examples:
            >>> human = db.with_species(" Homo sapiens ")  # doctest: +SKIP
            >>> human.species  # doctest: +SKIP
            'Homo sapiens'
        """

        normalized = str(species).strip()
        if not normalized:
            raise ValueError(
                "WikiPathways species must be non-empty after normalization"
            )
        if (
            self._publication_path is not None
            and self.snapshot.species is not None
            and normalized != self.snapshot.species
        ):
            raise CapabilityError(
                "WikiPathways publication is scoped to "
                f"{self.snapshot.species!r}; reopen an unscoped publication "
                "to query another species"
            )
        result = WikiPathwaysDatabase(
            snapshot=_WikiPathwaysSnapshot(
                files_gmt=self.snapshot.files_gmt,
                species=normalized,
            )
        )
        result._release_version = self._release_version
        result._publication_path = self._publication_path
        result._publication_identity = self._publication_identity
        return result

    def pathways(self) -> pl.LazyFrame:
        """Return one lazy metadata row per pathway in the current scope.

        Examples:
            >>> lf = db.pathways()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        return self._relation("pathway")

    def pathway_genes(self) -> pl.LazyFrame:
        """Return lazy pathway-to-Entrez membership in the current scope.

        Examples:
            >>> lf = db.pathway_genes()  # doctest: +SKIP
            >>> lf.select("term", "gene").collect()  # doctest: +SKIP
            shape: (..., 2)
        """

        return self._relation("term2gene")

    def pathway_names(self) -> pl.LazyFrame:
        """Return lazy display metadata for pathways in the current scope.

        Examples:
            >>> lf = db.pathway_names()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        return self._relation("term2name")

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> WikiPathwaysSelection:
        """Create a grouped selection from multiple Entrez Gene ID sets.

        Args:
            ids_by_group: Mapping from group label to input Entrez Gene IDs.

        Returns:
            A grouped selection that carries `group_id` through outputs.

        Raises:
            ValueError: If group IDs are invalid after normalization.

        Examples:
            Preserve comparison labels in the mapping table:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.select_groups({"case": ["2687"], "control": ["435"]})
            ...     .mappings().collect()
            ...     .select("group_id", "input_id", "wiki_pathways_id")
            ...     .to_dicts()
            ... )
            [{'group_id': 'case', 'input_id': '2687', 'wiki_pathways_id': 'WP100'}, {'group_id': 'control', 'input_id': '435', 'wiki_pathways_id': 'WP106'}]
        """
        self._assert_publication_current()
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return WikiPathwaysSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            _df_group_membership=grp_in_frames.df_group_membership,
        )

    def _eager_pathway(self) -> pl.DataFrame:
        """Materialize one metadata row per pathway in the current scope.

        ``gene_count`` counts distinct, non-empty Entrez Gene IDs in the GMT row.

        Examples:
            List pathway IDs in the compact human fixture:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ... ).with_species("Homo sapiens")
            >>> db.pathways().collect()["wiki_pathways_id"].to_list()
            ['WP100', 'WP106']
        """
        return self._lazy_frame("pathway").collect()

    def _eager_term2gene(self) -> pl.DataFrame:
        """Materialize distinct pathway-to-Entrez enrichment pairs.

        Examples:
            Inspect pathway-to-gene pairs from the compact human fixture:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ... ).with_species("Homo sapiens")
            >>> db.pathway_genes().collect().head(2).to_dicts()
            [{'wiki_pathways_id': 'WP100', 'gene_id': '2678'}, {'wiki_pathways_id': 'WP100', 'gene_id': '2687'}]
        """
        return self._lazy_frame("term2gene").collect()

    def _eager_term2name(self) -> pl.DataFrame:
        """Materialize one pathway display-metadata row per WikiPathways ID.

        Examples:
            Look up a display name for an enrichment term:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.pathway_names().collect()
            ...     .select("wiki_pathways_id", "pathway_name")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'wiki_pathways_id': 'WP100', 'pathway_name': 'Glutathione metabolism'}]
        """
        return self._lazy_frame("term2name").collect()

    def build_tidy(self) -> TidyDataset:
        """Build the lazy WikiPathways tidy dataset for the current scope.

        Returns:
            A `TidyDataset` with species-consistent ``pathway``, ``term2gene``,
            and ``term2name`` frames.

        Examples:
            Build the three declared WikiPathways frames:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> sorted(db.build_tidy().frames)
            ['pathway', 'term2gene', 'term2name']
        """
        self._assert_publication_current()
        frames = {
            "pathway": self._lazy_frame("pathway"),
            "term2gene": self._lazy_frame("term2gene"),
            "term2name": self._lazy_frame("term2name"),
        }
        release_version = self._release_version
        if release_version is None:
            raise RuntimeError("WikiPathways GMT release Version was not parsed")
        dataset_type = (
            _ReopenedWikiPathwaysTidyDataset
            if self._publication_path is not None
            else TidyDataset
        )
        return dataset_type(
            frames=frames,
            source=tuple(
                TidySource(
                    logical_name=f"pathway_gmt_{index:03d}",
                    path=file_gmt,
                    media_type=MEDIA_TYPE_GMT,
                )
                for index, file_gmt in enumerate(
                    self.snapshot.files_gmt,
                    start=1,
                )
            ),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="wikipathways-gmt-v1",
            build_id_prefix="wikipathways-gmt",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
            resource_name="wikipathways",
            release_version=release_version,
            release_version_source="official_metadata",
            scope=(
                json.dumps(
                    {"species": self.species},
                    separators=(",", ":"),
                    sort_keys=True,
                )
                if self.species is not None
                else None
            ),
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish pathway and membership relations as one DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(
            ...         Path(dir_out) / "wikipathways.duckdb"
            ...     )
            ...     result.tables
            ('pathway', 'pathway_gene')
        """
        if self._publication_path is not None:
            raise CapabilityError(
                "write_duckdb() requires a WikiPathways GMT source handle"
            )
        dataset = self.build_tidy()
        canonical = TidyDataset(
            frames=dataset.frames,
            source=dataset.source,
            resource_schema_version=dataset.resource_schema_version,
            source_schema_profile=dataset.source_schema_profile,
            build_id_prefix=dataset.build_id_prefix,
            assets=tuple(
                asset for asset in dataset.assets if asset.frame_name != "term2name"
            ),
            resource_name=dataset.resource_name,
            release_version=dataset.release_version,
            release_version_source=dataset.release_version_source,
            scope=dataset.scope,
        )
        return canonical.write_duckdb(
            path,
            table_names={"term2gene": "pathway_gene"},
            if_exists=if_exists,
        )

    def _lazy_frame(self, frame_name: str) -> pl.LazyFrame:
        """Return a cached lazy tidy frame in the current species scope.

        Args:
            frame_name: One of ``pathway``, ``term2gene``, or ``term2name``.

        Returns:
            The requested lazy frame. Parsing is deferred until the first frame
            request, then the frame mapping is reused by subsequent operations.

        Raises:
            KeyError: If ``frame_name`` is not a declared tidy frame.
            ValueError: If the GMT content is malformed when first parsed.

        Notes:
            Public extractors and `build_tidy()` own the supported frame
            contracts. Keep this cache entrypoint private so internal frame
            keys do not become a second API.
        """
        self._assert_publication_current()
        if self._frames is None:
            if self._publication_path is None:
                parsed = read_gmt_frames(self.snapshot.files_gmt)
                frames = parsed.frames
                self._release_version = parsed.release_version
            else:
                frames = self._read_publication_frames()
            if self.snapshot.species is not None:
                lf_pathway = _filter_species_frame(
                    frames["pathway"],
                    self.snapshot.species,
                )
                lf_pathway_ids = lf_pathway.select("wiki_pathways_id").unique()
                frames = {
                    "pathway": lf_pathway,
                    "term2gene": frames["term2gene"]
                    .join(
                        lf_pathway_ids,
                        on="wiki_pathways_id",
                        how="inner",
                    )
                    .sort("wiki_pathways_id", "gene_id"),
                    "term2name": _filter_species_frame(
                        frames["term2name"],
                        self.snapshot.species,
                    ),
                }
            self._frames = frames
        return self._frames[frame_name]

    def _relation(self, frame_name: str) -> pl.LazyFrame:
        schema = {
            "pathway": SCHEMA_PATHWAY,
            "term2gene": SCHEMA_TERM2GENE,
            "term2name": SCHEMA_TERM2NAME,
        }.get(frame_name)
        if schema is None:
            raise KeyError(f"Unknown WikiPathways relation: {frame_name}")
        if self._publication_path is None:
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_source_relation_batches(
                    self.snapshot.files_gmt,
                    species=self.snapshot.species,
                    frame_name=frame_name,
                    request=request,
                ),
            )
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_publication_relation_batches(
                self,
                frame_name=frame_name,
                request=request,
            ),
        )

    def _read_publication_frames(self) -> dict[str, pl.LazyFrame]:
        with self.connect() as connection:
            pathway = pl.read_database(  # pyright: ignore[reportUnknownMemberType]
                "SELECT * FROM pathway", connection
            )
            term2gene = pl.read_database(  # pyright: ignore[reportUnknownMemberType]
                "SELECT * FROM pathway_gene", connection
            )
        term2name = pathway.drop("gene_count")
        return {
            "pathway": pathway.lazy(),
            "term2gene": term2gene.lazy(),
            "term2name": term2name.lazy(),
        }

    def _assert_publication_identity(self) -> None:
        path = self._publication_path
        try:
            current_identity = None if path is None else _file_identity(path)
        except OSError:
            current_identity = None
        if current_identity != self._publication_identity:
            raise IntegrityError(
                "WikiPathways publication was replaced; reopen it with from_duckdb()"
            )

    def _assert_publication_current(self) -> None:
        if self._publication_path is not None:
            self._assert_publication_identity()


@dataclass(slots=True)
class WikiPathwaysSelection:
    """Selection handle for single and grouped WikiPathways queries.

    Selections are created by :meth:`WikiPathwaysDatabase.select_ids` or
    :meth:`WikiPathwaysDatabase.select_groups`. Single selections return tables keyed
    by `input_id`; grouped selections prepend `group_id`.

    Examples:
        Use a returned selection to materialize matched pathways:

        >>> db = WikiPathwaysDatabase.from_gmt(
        ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
        ... )
        >>> selection = db.select_ids(["2687"])
        >>> (
        ...     selection.mappings().collect()
        ...     .select("input_id", "wiki_pathways_id")
        ...     .to_dicts()
        ... )
        [{'input_id': '2687', 'wiki_pathways_id': 'WP100'}]
    """

    dataset: WikiPathwaysDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `group_id` through outputs.

        Examples:
            Inspect a grouped selection:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> selection = db.select_groups({"case": ["2687"]})
            >>> selection.is_grouped
            True
        """
        return self._df_groups is not None

    def mappings(self) -> pl.LazyFrame:
        """Return selected pathway mappings as a native lazy relation.

        Species scope is already applied to both the pathway metadata and its
        gene membership before this input join is planned.

        Examples:
            >>> lf = selection.mappings()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        snapshot = copy.copy(self)
        snapshot.dataset = copy.copy(self.dataset)
        schema = SCHEMA_GROUP_MAPPING if self.is_grouped else SCHEMA_MAPPING
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_selection_mapping_batches(
                snapshot,
                request=request,
            ),
        )

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return selected IDs absent from the current species scope.

        Examples:
            >>> lf = selection.unmatched_ids()  # doctest: +SKIP
            >>> lf.collect()  # doctest: +SKIP
            shape: (..., ...)
        """

        snapshot = copy.copy(self)
        schema = (
            SCHEMA_GROUP_INPUT_IDS
            if self._df_group_membership is not None
            else SCHEMA_UNMAPPED
        )
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_selection_unmatched_batches(
                snapshot,
                request=request,
            ),
        )


def _iter_selection_mapping_batches(
    selection: WikiPathwaysSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    batch_size = request.effective_batch_size
    database = selection.dataset
    if database._publication_path is None:  # pyright: ignore[reportPrivateUsage]
        requested = _requested_columns(
            request.columns,
            SCHEMA_GROUP_MAPPING if selection.is_grouped else SCHEMA_MAPPING,
        )
        input_ids = set(
            selection._df_input_ids.get_column("input_id").to_list()  # pyright: ignore[reportPrivateUsage]
        )
        membership = selection._df_group_membership  # pyright: ignore[reportPrivateUsage]
        rows: list[dict[str, object]] = []
        for pathway, gene_ids in iter_gmt_records(database.snapshot.files_gmt):
            if database.species is not None and pathway["species"] != database.species:
                continue
            for gene_id in gene_ids:
                if gene_id not in input_ids:
                    continue
                rows.append(
                    {
                        "input_id": gene_id,
                        "gene_id": gene_id,
                        "wiki_pathways_id": pathway["wiki_pathways_id"],
                        "pathway_name": pathway["pathway_name"],
                        "species": pathway["species"],
                        "url": pathway["url"],
                    }
                )
                if len(rows) >= batch_size:
                    yield from _finalize_source_mapping_batch(
                        rows,
                        membership=membership,
                        requested=requested,
                        batch_size=batch_size,
                    )
                    rows = []
        if rows:
            yield from _finalize_source_mapping_batch(
                rows,
                membership=membership,
                requested=requested,
                batch_size=batch_size,
            )
        return

    connection = database.connect()
    try:
        connection.execute(
            "CREATE TEMP TABLE _input_id(input_id VARCHAR NOT NULL PRIMARY KEY)"
        )
        input_rows = [
            (str(row["input_id"]),)
            for row in selection._df_input_ids.to_dicts()  # pyright: ignore[reportPrivateUsage]
        ]
        if input_rows:
            connection.executemany("INSERT INTO _input_id VALUES (?)", input_rows)
        query = """
            SELECT DISTINCT
                input.input_id AS input_id,
                input.input_id AS gene_id,
                gene.wiki_pathways_id AS wiki_pathways_id,
                pathway.pathway_name AS pathway_name,
                pathway.species AS species,
                pathway.url AS url
            FROM _input_id AS input
            JOIN pathway_gene AS gene
              ON gene.gene_id = input.input_id
            JOIN pathway
              ON pathway.wiki_pathways_id = gene.wiki_pathways_id
        """
        params: list[str] = []
        if database.species is not None:
            query += " WHERE pathway.species = ?"
            params.append(database.species)
        query += " ORDER BY input_id, wiki_pathways_id"
        requested = _requested_columns(
            request.columns,
            SCHEMA_GROUP_MAPPING if selection.is_grouped else SCHEMA_MAPPING,
        )
        result = connection.execute(query, params)
        reader = _publication_arrow_reader(result, batch_size)
        try:
            membership = selection._df_group_membership  # pyright: ignore[reportPrivateUsage]
            for record_batch in reader:
                frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
                if membership is not None:
                    frame = (
                        membership.join(frame, on="input_id", how="inner")
                        .select(list(SCHEMA_GROUP_MAPPING))
                        .unique()
                        .sort(list(SCHEMA_GROUP_MAPPING))
                    )
                if requested is not None:
                    frame = frame.select(requested)
                yield frame
        finally:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
    finally:
        connection.close()


def _iter_selection_unmatched_batches(
    selection: WikiPathwaysSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    database = selection.dataset
    input_ids = set(
        selection._df_input_ids.get_column("input_id").to_list()  # pyright: ignore[reportPrivateUsage]
    )
    if not input_ids:
        return
    mapped_ids: set[str] = set()
    if database._publication_path is None:  # pyright: ignore[reportPrivateUsage]
        for pathway, gene_ids in iter_gmt_records(database.snapshot.files_gmt):
            if database.species is not None and pathway["species"] != database.species:
                continue
            mapped_ids.update(input_ids.intersection(gene_ids))
    else:
        connection = database.connect()
        try:
            connection.execute(
                "CREATE TEMP TABLE _input_id(input_id VARCHAR NOT NULL PRIMARY KEY)"
            )
            connection.executemany(
                "INSERT INTO _input_id VALUES (?)",
                [(value,) for value in sorted(input_ids)],
            )
            query = (
                "SELECT DISTINCT input.input_id FROM _input_id AS input "
                "JOIN pathway_gene AS gene ON gene.gene_id = input.input_id "
                "JOIN pathway ON pathway.wiki_pathways_id = gene.wiki_pathways_id"
            )
            params: list[str] = []
            if database.species is not None:
                query += " WHERE pathway.species = ?"
                params.append(database.species)
            rows = connection.execute(query, params).fetchall()
            mapped_ids.update(str(row[0]) for row in rows)
        finally:
            connection.close()

    input_frame = selection._df_input_ids.filter(  # pyright: ignore[reportPrivateUsage]
        ~pl.col("input_id").is_in(mapped_ids)
    )
    if selection._df_group_membership is not None:  # pyright: ignore[reportPrivateUsage]
        input_frame = (
            selection._df_group_membership.join(  # pyright: ignore[reportPrivateUsage]
                input_frame,
                on="input_id",
                how="inner",
            )
            .select("group_id", "input_id")
            .sort("group_id", "input_id")
        )
    schema = (
        SCHEMA_GROUP_INPUT_IDS
        if selection._df_group_membership is not None  # pyright: ignore[reportPrivateUsage]
        else SCHEMA_UNMAPPED
    )
    requested = _requested_columns(request.columns, schema)
    if requested is not None:
        input_frame = input_frame.select(requested)
    for offset in range(0, input_frame.height, request.effective_batch_size):
        yield input_frame.slice(offset, request.effective_batch_size)


def _filter_species_frame(lf: pl.LazyFrame, species: str) -> pl.LazyFrame:
    if "species" not in lf.collect_schema().names():
        return lf
    return lf.filter(pl.col("species") == species)


def _finalize_source_mapping_batch(
    rows: list[dict[str, object]],
    *,
    membership: pl.DataFrame | None,
    requested: list[str] | None,
    batch_size: int,
) -> Iterator[pl.DataFrame]:
    frame = pl.DataFrame(rows, schema=SCHEMA_MAPPING)
    if membership is not None:
        frame = (
            membership.join(frame, on="input_id", how="inner")
            .select(list(SCHEMA_GROUP_MAPPING))
            .unique()
            .sort(list(SCHEMA_GROUP_MAPPING))
        )
    else:
        frame = frame.unique().sort(list(SCHEMA_MAPPING))
    if requested is not None:
        frame = frame.select(requested)
    for offset in range(0, frame.height, batch_size):
        yield frame.slice(offset, batch_size)


def _iter_source_relation_batches(
    files_gmt: tuple[Path, ...],
    *,
    species: str | None,
    frame_name: str,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    requested = _requested_columns(
        request.columns,
        {
            "pathway": SCHEMA_PATHWAY,
            "term2gene": SCHEMA_TERM2GENE,
            "term2name": SCHEMA_TERM2NAME,
        }[frame_name],
    )
    schema = {
        "pathway": SCHEMA_PATHWAY,
        "term2gene": SCHEMA_TERM2GENE,
        "term2name": SCHEMA_TERM2NAME,
    }[frame_name]
    rows: list[dict[str, object]] = []
    batch_size = request.effective_batch_size
    for pathway, gene_ids in iter_gmt_records(files_gmt):
        if species is not None and pathway["species"] != species:
            continue
        if frame_name == "pathway":
            rows.append(dict(pathway))
        elif frame_name == "term2gene":
            for gene_id in gene_ids:
                rows.append(
                    {
                        "wiki_pathways_id": pathway["wiki_pathways_id"],
                        "gene_id": gene_id,
                    }
                )
                if len(rows) >= batch_size:
                    frame = pl.DataFrame(rows, schema=schema)
                    if requested is not None:
                        frame = frame.select(requested)
                    yield frame
                    rows = []
        else:
            rows.append(
                {
                    "wiki_pathways_id": pathway["wiki_pathways_id"],
                    "pathway_name": pathway["pathway_name"],
                    "species": pathway["species"],
                    "collection": pathway["collection"],
                    "version": pathway["version"],
                    "url": pathway["url"],
                }
            )
        if len(rows) >= batch_size:
            frame = pl.DataFrame(rows, schema=schema)
            if requested is not None:
                frame = frame.select(requested)
            yield frame
            rows = []
    if rows:
        frame = pl.DataFrame(rows, schema=schema)
        if requested is not None:
            frame = frame.select(requested)
        yield frame


def _iter_publication_relation_batches(
    database: WikiPathwaysDatabase,
    *,
    frame_name: str,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    batch_size = request.effective_batch_size
    species = database.species
    if frame_name == "pathway":
        columns = _requested_columns(request.columns, SCHEMA_PATHWAY)
        query = "SELECT " + _select_sql(columns, SCHEMA_PATHWAY) + " FROM pathway"
        params: list[str] = []
        if species is not None:
            query += " WHERE species = ?"
            params.append(species)
    elif frame_name == "term2gene":
        columns = _requested_columns(request.columns, SCHEMA_TERM2GENE)
        query = """
            SELECT gene.wiki_pathways_id, gene.gene_id
            FROM pathway_gene AS gene
        """
        params = []
        if species is not None:
            query += " JOIN pathway AS pathway USING (wiki_pathways_id)"
            query += " WHERE pathway.species = ?"
            params.append(species)
        query += " ORDER BY gene.wiki_pathways_id, gene.gene_id"
    elif frame_name == "term2name":
        columns = _requested_columns(request.columns, SCHEMA_TERM2NAME)
        query = """
            SELECT
                wiki_pathways_id, pathway_name, species,
                collection, version, url
            FROM pathway
        """
        params = []
        if species is not None:
            query += " WHERE species = ?"
            params.append(species)
        query += " ORDER BY wiki_pathways_id"
    else:
        raise KeyError(f"Unknown WikiPathways relation: {frame_name}")

    connection = database.connect()
    try:
        if columns is not None:
            query = _project_query(query, columns)
        result = connection.execute(query, params)
        reader = _publication_arrow_reader(result, batch_size)
        try:
            for record_batch in reader:
                frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
                yield frame
        finally:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
    finally:
        connection.close()


def _requested_columns(
    columns: tuple[str, ...] | None,
    schema: SchemaDict,
) -> list[str] | None:
    if columns is None:
        return None
    selected = [name for name in columns if name in schema]
    return selected or None


def _select_sql(columns: list[str] | None, schema: SchemaDict) -> str:
    selected = list(schema) if columns is None else columns
    return ", ".join(f'"{name}"' for name in selected)


def _project_query(query: str, columns: list[str]) -> str:
    projection = ", ".join(f'"{name}"' for name in columns)
    return f"SELECT {projection} FROM ({query}) AS _bioextract_relation"


def _publication_arrow_reader(result: Any, batch_size: int) -> Any:
    to_arrow_reader = getattr(result, "to_arrow_reader", None)
    if to_arrow_reader is not None:
        return to_arrow_reader(batch_size)
    fetch_record_batch = result.fetch_record_batch
    return fetch_record_batch(rows_per_batch=batch_size)


_WIKIPATHWAYS_TABLE_CONTRACTS = {
    "pathway": (
        "canonical",
        (
            ("wiki_pathways_id", "VARCHAR"),
            ("pathway_name", "VARCHAR"),
            ("species", "VARCHAR"),
            ("collection", "VARCHAR"),
            ("version", "VARCHAR"),
            ("url", "VARCHAR"),
            ("gene_count", "BIGINT"),
        ),
    ),
    "pathway_gene": (
        "derived",
        (("wiki_pathways_id", "VARCHAR"), ("gene_id", "VARCHAR")),
    ),
}


def _validate_wikipathways_publication(path: Path) -> tuple[str, str | None]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with duckdb.connect(str(path), read_only=True) as connection:
            metadata_rows = connection.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
            metadata = {str(row[0]): str(row[1]) for row in metadata_rows}
            if len(metadata) != len(metadata_rows):
                raise ValueError("WikiPathways publication has duplicate metadata keys")
            if (
                metadata.get("bioextract.metadata_schema_version")
                != METADATA_SCHEMA_VERSION
            ):
                raise ValueError("Unsupported WikiPathways metadata schema version")
            validate_duckdb_metadata_v2(connection, metadata)
            if metadata.get("bioextract.resource_name") != "wikipathways":
                raise ValueError(
                    "DuckDB file is not a bioextract WikiPathways publication"
                )
            if metadata.get("bioextract.source_schema_profile") != (
                "wikipathways-gmt-v1"
            ):
                raise ValueError("Unsupported WikiPathways source schema profile")
            if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported WikiPathways resource schema version")
            declared_species: str | None = None
            scope_value = metadata.get("bioextract.scope")
            if scope_value is not None:
                scope_obj = json.loads(scope_value)
                if not isinstance(scope_obj, dict):
                    raise ValueError("WikiPathways publication scope is unsupported")
                scope = cast(dict[str, object], scope_obj)
                if set(scope) != {"species"}:
                    raise ValueError("WikiPathways publication scope is unsupported")
                declared_species = str(scope["species"]).strip()
                if not declared_species:
                    raise ValueError("WikiPathways publication species scope is empty")
            release_version = metadata.get("bioextract.release_version")
            if (
                release_version is None
                or metadata.get("bioextract.release_version_source")
                != "official_metadata"
            ):
                raise ValueError(
                    "WikiPathways publication release identity is unsupported"
                )

            source_rows = connection.execute(
                "SELECT logical_name, bytes, media_type FROM _bioextract.source_file"
            ).fetchall()
            expected_roles = {
                f"pathway_gmt_{index:03d}" for index in range(1, len(source_rows) + 1)
            }
            if (
                not source_rows
                or {str(row[0]) for row in source_rows} != expected_roles
                or any(
                    (row[1] is not None and int(row[1]) < 0)
                    or str(row[2]) != MEDIA_TYPE_GMT
                    for row in source_rows
                )
            ):
                raise ValueError(
                    "WikiPathways source capability inventory is unsupported"
                )

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
            } | {("main", name, "BASE TABLE") for name in _WIKIPATHWAYS_TABLE_CONTRACTS}
            if relations != expected_relations:
                raise ValueError(
                    "WikiPathways physical table/view inventory is unsupported"
                )

            info_rows = connection.execute(
                "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
            ).fetchall()
            recorded = {str(row[0]): (str(row[1]), int(row[2])) for row in info_rows}
            if len(recorded) != len(info_rows) or set(recorded) != set(
                _WIKIPATHWAYS_TABLE_CONTRACTS
            ):
                raise ValueError("WikiPathways table inventory does not match metadata")
            for table_name, (role, row_count) in recorded.items():
                expected_role, expected_schema = _WIKIPATHWAYS_TABLE_CONTRACTS[
                    table_name
                ]
                if role != expected_role or row_count < 0:
                    raise ValueError(
                        "WikiPathways table capability inventory is unsupported"
                    )
                actual_schema = tuple(
                    (str(row[1]), str(row[2]))
                    for row in connection.execute(
                        f'PRAGMA table_info("{table_name}")'
                    ).fetchall()
                )
                if actual_schema != expected_schema:
                    raise ValueError(
                        f"WikiPathways table schema is unsupported: {table_name}"
                    )
            column_mappings = {
                tuple(str(value) for value in row)
                for row in connection.execute(
                    "SELECT table_name, source_column, output_column, reason "
                    "FROM _bioextract.column_mapping"
                ).fetchall()
            }
            if column_mappings:
                raise ValueError(
                    "WikiPathways column provenance inventory is unsupported"
                )
    except duckdb.Error as error:
        raise ValueError(
            f"Cannot open WikiPathways DuckDB publication: {path}"
        ) from error
    return release_version, declared_species


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
