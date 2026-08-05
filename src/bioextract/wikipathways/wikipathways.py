from __future__ import annotations

import os
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import polars as pl

from bioextract._publication import (
    BIOEXTRACT_RELATIONS,
    DuckDBWriteResult,
    validate_duckdb_metadata_v1,
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
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
)
from .util import (
    extract_mapping_frame,
    extract_unmatched_ids_frame,
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
        include_source_hashes: bool = False,
    ) -> DuckDBWriteResult:
        del path, table_names, if_exists, source_columns, include_source_hashes
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
        ...     species="Homo sapiens",
        ... )
        >>> selection = db.select_ids(["2687", "MISSING"])
        >>> (
        ...     selection.extract_mapping()
        ...     .select("input_id", "wiki_pathways_id")
        ...     .to_dicts()
        ... )
        [{'input_id': '2687', 'wiki_pathways_id': 'WP100'}]
        >>> selection.extract_unmatched_ids().to_dicts()
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
        species: str | None = None,
        glob: bool = True,
    ) -> WikiPathwaysDatabase:
        """Create a dataset handle from one or more local WikiPathways GMT files.

        Args:
            source: One local path, a sequence of local paths, or glob
                expressions. Every resolved physical file must be unique.
            species: Optional species display name used as an exact metadata
                filter after parsing.
            glob: Expand each source entry as a glob expression. If false,
                treat every entry literally.

        Returns:
            A dataset handle that can build pathway, term2gene, and term2name
            frames.

        Raises:
            FileNotFoundError: If a literal file does not exist or a glob
                expression has no matches.
            ValueError: At construction, if the source is empty, resolves to a
                directory, repeats a physical file, or the species string is
                empty after normalization. GMT content parsing is deferred;
                malformed content, inconsistent Collection or Version values,
                and duplicate WikiPathways IDs raise when a frame is first
                extracted, built, or written.

        Examples:
            Open all GMT files in one local snapshot and retain human rows:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/*.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> (
            ...     db.extract_pathway()
            ...     .select("wiki_pathways_id", "pathway_name")
            ...     .head(1)
            ...     .to_dicts()
            ... )
            [{'wiki_pathways_id': 'WP100', 'pathway_name': 'Glutathione metabolism'}]
        """
        files_gmt = resolve_gmt_sources(source, glob=glob)
        species_normalized = None if species is None else str(species).strip()
        if species is not None and not species_normalized:
            raise ValueError(
                "WikiPathways species must be non-empty after normalization"
            )
        return cls(
            snapshot=_WikiPathwaysSnapshot(
                files_gmt=files_gmt,
                species=species_normalized,
            ),
        )

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> WikiPathwaysDatabase:
        """Open a validated WikiPathways publication for domain and SQL access.

        Validation reads bounded metadata and catalog schemas only. The handle
        is pinned to the exact file that passed validation.

        Args:
            path: A bioextract WikiPathways metadata-v1 DuckDB publication.

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
            >>> db.extract_pathway().height > 0  # doctest: +SKIP
            True
        """
        publication_path = Path(path).absolute()
        identity_before = _file_identity(publication_path)
        try:
            release_version = _validate_wikipathways_publication(publication_path)
            identity_after = _file_identity(publication_path)
        except (KeyError, TypeError, ValueError) as error:
            raise IntegrityError(str(error)) from error
        except OSError as error:
            raise IntegrityError(
                "WikiPathways publication changed during validation"
            ) from error
        if identity_after != identity_before:
            raise IntegrityError("WikiPathways publication changed during validation")
        result = cls(snapshot=_WikiPathwaysSnapshot(files_gmt=()))
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
            ...     .extract_mapping()
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
            ...     .extract_mapping()
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

    def extract_pathway(self) -> pl.DataFrame:
        """Extract one metadata row per pathway in the current species scope.

        ``gene_count`` counts distinct, non-empty Entrez Gene IDs in the GMT row.

        Examples:
            List pathway IDs in the compact human fixture:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_pathway()["wiki_pathways_id"].to_list()
            ['WP100', 'WP106']
        """
        return self._lazy_frame("pathway").collect()

    def extract_term2gene(self) -> pl.DataFrame:
        """Extract distinct WikiPathways-pathway-to-Entrez enrichment pairs.

        Examples:
            Inspect pathway-to-gene pairs from the compact human fixture:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> db.extract_term2gene().head(2).to_dicts()
            [{'wiki_pathways_id': 'WP100', 'gene_id': '2678'}, {'wiki_pathways_id': 'WP100', 'gene_id': '2687'}]
        """
        return self._lazy_frame("term2gene").collect()

    def extract_term2name(self) -> pl.DataFrame:
        """Extract one pathway display-metadata row per WikiPathways ID.

        Examples:
            Look up a display name for an enrichment term:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> (
            ...     db.extract_term2name()
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
        ...     selection.extract_mapping()
        ...     .select("input_id", "wiki_pathways_id")
        ...     .to_dicts()
        ... )
        [{'input_id': '2687', 'wiki_pathways_id': 'WP100'}]
    """

    dataset: WikiPathwaysDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

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

    def extract_mapping(self) -> pl.DataFrame:
        """Extract every pathway mapping matched by the selected Entrez IDs.

        Examples:
            Materialize the pathway mapped to Entrez ID 2687:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt",
            ...     species="Homo sapiens",
            ... )
            >>> selection = db.select_ids(["2687"])
            >>> selection.extract_mapping()["wiki_pathways_id"].to_list()
            ['WP100']
        """
        self.dataset._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        if self._df_mapping is None:
            self._df_mapping = extract_mapping_frame(
                self.dataset._lazy_frame("pathway"),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                self.dataset._lazy_frame("term2gene"),  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
                self._df_input_ids,
                df_group_membership=self._df_group_membership,
            ).collect()
        return self._df_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no WikiPathways mapping.

        Grouped selections report an ID as unmapped independently within each
        group and include ``group_id`` in the result.

        Examples:
            Retain an Entrez ID absent from the snapshot:

            >>> db = WikiPathwaysDatabase.from_gmt(
            ...     "data/wikipathways/wikipathways-Homo_sapiens.gmt"
            ... )
            >>> selection = db.select_ids(["2687", "MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'input_id': 'MISSING'}]
        """
        self.dataset._assert_publication_current()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                df_group_membership=self._df_group_membership,
            )
        return self._df_unmapped


def _filter_species_frame(lf: pl.LazyFrame, species: str) -> pl.LazyFrame:
    if "species" not in lf.collect_schema().names():
        return lf
    return lf.filter(pl.col("species") == species)


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


def _validate_wikipathways_publication(path: Path) -> str:
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
            if metadata.get("bioextract.metadata_schema_version") != "1":
                raise ValueError("Unsupported WikiPathways metadata schema version")
            validate_duckdb_metadata_v1(connection, metadata)
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
            if "bioextract.scope" in metadata:
                raise ValueError("WikiPathways publication scope is unsupported")
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
                    int(row[1]) < 0 or str(row[2]) != MEDIA_TYPE_GMT
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
    return release_version


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
