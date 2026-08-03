from __future__ import annotations

import json
import os
from collections.abc import Iterable, Mapping
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
from bioextract._tidy import (
    TidyAsset,
    TidyDataset,
    TidySource,
)

from ._pfam import build_pfam_tidy_dataset, read_interpro_release_version
from .constant import (
    ASSET_SPECS,
    MEDIA_TYPE_TSV_GZIP,
    MEDIA_TYPE_XML_GZIP,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    SCHEMA_VERSION,
    InterProNamespace,
)
from .util import (
    extract_mapping_frame,
    extract_unmatched_ids_frame,
    read_interpro_xml_frames,
    read_mapping_frame,
    scan_mapping_frame,
    select_mapping_frame,
    validate_mapping_xml_relationships,
)

__all__ = [
    "InterProDatabase",
]


@dataclass(frozen=True, slots=True)
class _InterProSnapshot:
    file_protein2ipr: Path | None = None
    file_interpro_xml: Path | None = None


@dataclass(slots=True)
class InterProDatabase:
    """Access one local InterPro protein-mapping snapshot.

    `protein2ipr.dat.gz` supplies the canonical mapping. Optional InterPro XML
    enriches mapping rows with entry type and member-database metadata; without
    XML those fields remain null. The compact Pfam tidy configuration validates
    exact `InterProId + PfamId` relationships against the supplied XML.

    Examples:
        Read the first domain annotation from a versioned snapshot:

        >>> db = InterProDatabase.from_mapping_files(
        ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
        ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
        ... )
        >>> db.extract_mapping().select(
        ...     "UniProtId", "InterProId", "MemberDb"
        ... ).head(1).rows()
        [('P12345', 'IPR000001', 'PFAM')]
    """

    snapshot: _InterProSnapshot
    _df_mapping: pl.DataFrame | None = field(default=None, init=False, repr=False)
    _frames_xml: dict[str, pl.DataFrame] | None = field(
        default=None, init=False, repr=False
    )
    _publication_path: Path | None = field(default=None, init=False, repr=False)
    _capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset[str](), init=False
    )
    release_version: str | None = field(default=None, init=False)

    @classmethod
    def from_mapping_files(
        cls,
        *,
        protein_to_interpro: os.PathLike[str] | str,
        interpro_xml: os.PathLike[str] | str | None = None,
    ) -> InterProDatabase:
        """Create a handle from local InterPro mapping files.

        Args:
            protein_to_interpro: Path to the required `protein2ipr.dat.gz` source.
            interpro_xml: Optional `interpro.xml.gz` source for mapping
                enrichment and required Pfam metadata.

        Returns:
            A lightweight handle that defers parsing until data is requested.

        Raises:
            FileNotFoundError: If a supplied source file does not exist.

        Examples:
            Use XML to recover entry and member-database metadata:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> db.extract_mapping().select(
            ...     "InterProId", "InterProType", "MemberDb"
            ... ).head(1).rows()
            [('IPR000001', 'Domain', 'PFAM')]
        """
        file_protein2ipr = _validate_file(
            protein_to_interpro,
            label="InterPro protein2ipr file",
        )
        file_interpro_xml = interpro_xml
        if file_interpro_xml is not None:
            file_interpro_xml = _validate_file(
                file_interpro_xml,
                label="InterPro XML file",
            )
        return cls(
            snapshot=_InterProSnapshot(
                file_protein2ipr=file_protein2ipr,
                file_interpro_xml=file_interpro_xml,
            ),
        )

    @classmethod
    def from_duckdb(cls, path: os.PathLike[str] | str) -> InterProDatabase:
        """Open a validated InterPro publication for domain and SQL access.

        Examples:
            >>> db = InterProDatabase.from_duckdb(  # doctest: +SKIP
            ...     "tidy/data.duckdb"
            ... )
            >>> db.select_ids(["P12345"]).extract_mapping().height  # doctest: +SKIP
            1
        """
        publication_path = Path(path)
        capabilities, release_version = _validate_interpro_publication(publication_path)
        result = cls(snapshot=_InterProSnapshot())
        result._publication_path = publication_path
        result._capabilities = capabilities
        result.release_version = release_version
        return result

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a fresh caller-owned read-only publication connection.

        Examples:
            >>> db = InterProDatabase.from_duckdb("tidy/data.duckdb")  # doctest: +SKIP
            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.sql("SELECT count(*) FROM mapping").fetchone()
            ...     count[0] >= 0
            True
        """
        if self._publication_path is None:
            raise ValueError("connect() requires InterProDatabase.from_duckdb()")
        return duckdb.connect(str(self._publication_path), read_only=True)

    def select_ids(
        self,
        ids: Iterable[str],
    ) -> InterProSelection:
        """Create a selection for one normalized UniProt ID set.

        Args:
            ids: UniProt accessions. Empty values are discarded and duplicate
                normalized IDs collapse to one input row.
        Returns:
            A selection handle whose mapping output includes `InputId` and
            `InputNamespace` provenance columns.

        Examples:
            Normalize a UniProt entry label before selecting its mapping rows:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(["sp|P12345|TEST_HUMAN"])
            >>> selection.extract_mapping().select(
            ...     "InputId", "InterProId", "MemberDbId"
            ... ).rows()
            [('P12345', 'IPR000001', 'PF00051'), ('P12345', 'IPR000001', 'SM00130')]
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return InterProSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            _df_group_membership=None,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> InterProSelection:
        """Create a selection that preserves caller-defined groups.

        Args:
            ids_by_group: Mapping from unique, non-empty group IDs to UniProt
                accessions.
        Returns:
            A selection handle whose mapping and unmapped outputs retain
            `GroupId`.

        Raises:
            ValueError: If group IDs are invalid.

        Examples:
            Preserve comparison labels on every selected mapping row:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_groups(
            ...     {"up": ["P12345"], "down": ["Q9Y243"]},
            ... )
            >>> selection.extract_mapping().select(
            ...     "GroupId", "UniProtId"
            ... ).unique().sort("GroupId").rows()
            [('down', 'Q9Y243'), ('up', 'P12345')]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return InterProSelection(
            dataset=self,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_groups=grp_in_frames.df_groups,
            _df_group_membership=grp_in_frames.df_group_membership,
        )

    def extract_mapping(self) -> pl.DataFrame:
        """Extract the full row-level UniProt-to-InterPro mapping.

        Returns:
            A cached frame that preserves member-database rows and positional
            coordinates. `InterProType` and `MemberDb` remain null when no XML
            source was configured.

        Examples:
            Retain each member-database match and its coordinates:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> db.extract_mapping().select(
            ...     "UniProtId", "MemberDbId", "Start", "End"
            ... ).head(2).rows()
            [('P12345', 'PF00051', 10, 80), ('P12345', 'SM00130', 12, 76)]
        """
        if self._df_mapping is None:
            if self._publication_path is not None:
                with self.connect() as connection:
                    self._df_mapping = pl.read_database(  # pyright: ignore[reportUnknownMemberType]  # Polars-DuckDB boundary
                        "SELECT * FROM mapping", connection
                    ).rename(
                        {name: source for source, name in _MAPPING_COLUMNS.items()}
                    )
            else:
                self._df_mapping = read_mapping_frame(
                    self._require_mapping_source(action="extract mappings"),
                    df_interpro_entry=self.xml_frame("entry"),
                    df_interpro_member=self.xml_frame("member"),
                )
        return self._df_mapping

    def build_tidy(
        self,
    ) -> TidyDataset:
        """Build all lazy relations supported by the configured source files.

        Returns:
            A mapping dataset, plus compact Pfam relations when XML is present.

        Raises:
            ValueError: If snapshot metadata and raw mappings fail validation.

        Notes:
            Mapping accepts missing XML and leaves enrichment fields null.
            Pfam relations are added when XML is available. Paths never supply
            release identity; official XML metadata does.

        Examples:
            The XML below includes an INTERPRO release record whose official
            version is published as the release identity.

            Build the default mapping plan:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> sorted(db.build_tidy().frames)
            ['mapping', 'protein_term', 'term', 'term_xref']
        """
        file_protein2ipr = self._require_mapping_source(action="build a publication")

        release_version = (
            read_interpro_release_version(self.snapshot.file_interpro_xml)
            if self.snapshot.file_interpro_xml is not None
            else None
        )
        if self.snapshot.file_interpro_xml is not None:
            validate_mapping_xml_relationships(
                file_protein2ipr,
                df_interpro_entry=self.xml_frame("entry"),
                df_interpro_member=self.xml_frame("member"),
            )
        dataset = TidyDataset(
            frames={
                "mapping": scan_mapping_frame(
                    file_protein2ipr,
                    df_interpro_entry=self.xml_frame("entry"),
                    df_interpro_member=self.xml_frame("member"),
                )
            },
            source=self._tidy_sources(),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile=(
                "interpro-protein2ipr-xml-v1"
                if self.snapshot.file_interpro_xml is not None
                else "interpro-protein2ipr-v1"
            ),
            build_id_prefix=f"interpro-{file_protein2ipr.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
            resource_name="interpro",
            release_version=release_version,
            release_version_source=(
                "official_metadata" if release_version is not None else None
            ),
            extra_metadata={
                "bioextract.capabilities": (
                    "mapping,pfam"
                    if self.snapshot.file_interpro_xml is not None
                    else "mapping"
                )
            },
        )
        if self.snapshot.file_interpro_xml is None:
            return dataset
        pfam = build_pfam_tidy_dataset(
            file_protein2ipr=file_protein2ipr,
            file_interpro_xml=self.snapshot.file_interpro_xml,
        )
        dataset.frames = {**dataset.frames, **pfam.frames}
        dataset.assets = (
            *dataset.assets,
            *(
                TidyAsset(asset.path, "compact", asset.frame_name)
                for asset in pfam.assets
            ),
        )
        return dataset

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
        include_source_hashes: bool = False,
    ) -> DuckDBWriteResult:
        """Publish every relation available from this source handle as DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(
            ...         Path(dir_out) / "interpro.duckdb"
            ...     )
            ...     result.tables
            ('mapping', 'protein_term', 'term', 'term_xref')
        """
        return self.build_tidy().write_duckdb(
            path,
            if_exists=if_exists,
            include_source_hashes=include_source_hashes,
        )

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

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> db.xml_frame("entry").rows()
            [('IPR000001', 'Domain'), ('IPR000002', 'Homologous_superfamily')]
        """
        if self._frames_xml is None:
            self._frames_xml = read_interpro_xml_frames(self.snapshot.file_interpro_xml)
        return self._frames_xml[frame_name]

    def _tidy_sources(self) -> tuple[TidySource, ...]:
        sources = [
            TidySource(
                logical_name="protein_to_interpro",
                path=self._require_mapping_source(
                    action="describe publication sources"
                ),
                media_type=MEDIA_TYPE_TSV_GZIP,
            )
        ]
        if self.snapshot.file_interpro_xml is not None:
            sources.append(
                TidySource(
                    logical_name="interpro_xml",
                    path=self.snapshot.file_interpro_xml,
                    media_type=MEDIA_TYPE_XML_GZIP,
                )
            )
        return tuple(sources)

    def _require_interpro_xml(self, *, action: str) -> Path:
        if self.snapshot.file_interpro_xml is None:
            raise ValueError(f"InterPro XML file is required to {action}")
        return self.snapshot.file_interpro_xml

    def _require_mapping_source(self, *, action: str) -> Path:
        if self.snapshot.file_protein2ipr is None:
            raise ValueError(f"InterPro mapping source is required to {action}")
        return self.snapshot.file_protein2ipr


@dataclass(slots=True)
class InterProSelection:
    """Materialize one single or grouped InterPro mapping query.

    Instances are created by `InterProDatabase.select_ids()` or
    `InterProDatabase.select_groups()`. Mapping and unmapped frames are cached
    independently after first extraction.

    Examples:
        Inspect the InterPro IDs retained by a selection:

        >>> db = InterProDatabase.from_mapping_files(
        ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
        ... )
        >>> selection = db.select_ids(["P12345"])
        >>> selection.extract_mapping().get_column("InterProId").unique().to_list()
        ['IPR000001']
    """

    dataset: InterProDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    namespace: InterProNamespace = field(default="uniprot", init=False)
    _df_mapping: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> db.select_groups(
            ...     {"up": ["P12345"]}
            ... ).is_grouped
            True
        """
        return self._df_groups is not None

    def extract_mapping(self) -> pl.DataFrame:
        """Extract mapping rows for the normalized selection.

        Returns:
            A frame prefixed by `InputId` and `InputNamespace`; grouped selections
            additionally begin with `GroupId`. Member and positional rows stay
            expanded.

        Examples:
            Extract the normalized input alongside its member-database rows:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(["P12345"])
            >>> selection.extract_mapping().select(
            ...     "InputId", "MemberDbId"
            ... ).rows()
            [('P12345', 'PF00051'), ('P12345', 'SM00130')]
        """
        if self._df_mapping is None:
            if self.dataset._publication_path is not None:  # pyright: ignore[reportPrivateUsage]
                selected = extract_mapping_frame(
                    self.dataset.extract_mapping(),
                    self._df_input_ids,
                    namespace=self.namespace,
                )
                if self._df_group_membership is None:
                    self._df_mapping = selected
                else:
                    columns = ["GroupId", *selected.columns]
                    self._df_mapping = (
                        self._df_group_membership.join(
                            selected, on="InputId", how="inner"
                        )
                        .select(columns)
                        .unique()
                        .sort(columns)
                    )
            else:
                self._df_mapping = select_mapping_frame(
                    self.dataset._require_mapping_source(  # pyright: ignore[reportPrivateUsage]
                        action="select mappings"
                    ),
                    self._df_input_ids,
                    df_group_membership=self._df_group_membership,
                    namespace=self.namespace,
                    df_interpro_entry=self.dataset.xml_frame("entry"),
                    df_interpro_member=self.dataset.xml_frame("member"),
                )
        return self._df_mapping

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs with no InterPro mapping row.

        Returns:
            `InputId` for a single selection, or `GroupId, InputId` for a
            grouped selection.

        Examples:
            Report an accession absent from the local snapshot:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> selection = db.select_ids(["MISSING"])
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            self._df_unmapped = extract_unmatched_ids_frame(
                self._df_input_ids,
                self.extract_mapping(),
                df_group_membership=self._df_group_membership,
            )
        return self._df_unmapped


def _validate_file(
    file_path: os.PathLike[str] | str,
    *,
    label: str,
) -> Path:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} not found: {file_path}")
    return file_path


_MAPPING_COLUMNS = {
    "UniProtId": "uniprot_id",
    "InterProId": "interpro_id",
    "InterProName": "interpro_name",
    "InterProType": "interpro_type",
    "MemberDb": "member_db",
    "MemberDbId": "member_db_id",
    "Start": "start",
    "End": "end",
}
_TABLE_CONTRACTS = {
    "mapping": (
        "canonical",
        [
            (name, "BIGINT" if name in {"start", "end"} else "VARCHAR")
            for name in _MAPPING_COLUMNS.values()
        ],
    ),
    "protein_term": (
        "compact",
        [("uniprot_id", "VARCHAR"), ("pfam_id", "VARCHAR")],
    ),
    "term": (
        "compact",
        [("pfam_id", "VARCHAR"), ("pfam_name", "VARCHAR")],
    ),
    "term_xref": (
        "compact",
        [
            ("pfam_id", "VARCHAR"),
            ("interpro_id", "VARCHAR"),
            ("interpro_name", "VARCHAR"),
            ("interpro_type", "VARCHAR"),
        ],
    ),
}


def _validate_interpro_publication(path: Path) -> tuple[frozenset[str], str | None]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with duckdb.connect(str(path), read_only=True) as connection:
            metadata = dict(
                connection.execute(
                    "SELECT key, value FROM _bioextract.metadata"
                ).fetchall()
            )
            if metadata.get("bioextract.metadata_schema_version") != "1":
                raise ValueError("Unsupported InterPro metadata schema version")
            validate_duckdb_metadata_v1(connection, metadata)
            if metadata.get("bioextract.resource_name") != "interpro":
                raise ValueError("DuckDB file is not a bioextract InterPro publication")
            if metadata.get("bioextract.resource_schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported InterPro resource schema version")

            capabilities = frozenset(
                value.strip()
                for value in metadata.get("bioextract.capabilities", "").split(",")
                if value.strip()
            )
            profile_value = metadata.get("bioextract.source_schema_profile")
            if not isinstance(profile_value, str):
                raise ValueError(
                    "InterPro publication is missing source schema profile"
                )
            profile_contracts = {
                "interpro-protein2ipr-v1": (
                    frozenset({"mapping"}),
                    {"mapping"},
                    {"protein_to_interpro"},
                ),
                "interpro-protein2ipr-xml-v1": (
                    frozenset({"mapping", "pfam"}),
                    set(_TABLE_CONTRACTS),
                    {"protein_to_interpro", "interpro_xml"},
                ),
            }
            expected = profile_contracts.get(profile_value)
            if expected is None:
                raise ValueError("Unsupported InterPro source schema profile")
            expected_capabilities, expected_tables, expected_sources = expected
            if capabilities != expected_capabilities:
                raise ValueError("InterPro capability inventory is unsupported")

            source_rows = connection.execute(
                "SELECT logical_name, display_path, bytes, media_type, sha256 "
                "FROM _bioextract.source_file ORDER BY logical_name"
            ).fetchall()
            source_roles = {str(row[0]) for row in source_rows}
            if source_roles != expected_sources:
                raise ValueError("InterPro source role inventory is unsupported")
            embedded_sources = json.loads(metadata["bioextract.sources"])
            table_sources = [
                {
                    "logical_name": str(row[0]),
                    "path": str(row[1]),
                    "bytes": int(row[2]),
                    "media_type": str(row[3]),
                    **({"sha256": str(row[4])} if row[4] is not None else {}),
                }
                for row in source_rows
            ]
            if (
                sorted(embedded_sources, key=lambda item: item["logical_name"])
                != table_sources
            ):
                raise ValueError("InterPro embedded source inventory is unsupported")
            expected_media_types = {
                "protein_to_interpro": MEDIA_TYPE_TSV_GZIP,
                "interpro_xml": MEDIA_TYPE_XML_GZIP,
            }
            if any(row[3] != expected_media_types[str(row[0])] for row in source_rows):
                raise ValueError("InterPro source media-type inventory is unsupported")

            provenance_tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT table_name FROM information_schema.tables "
                    "WHERE table_schema='_bioextract' AND table_type='BASE TABLE'"
                ).fetchall()
            }
            if provenance_tables != BIOEXTRACT_RELATIONS:
                raise ValueError("InterPro provenance table inventory is unsupported")
            physical_tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT table_name FROM information_schema.tables "
                    "WHERE table_schema='main' AND table_type='BASE TABLE'"
                ).fetchall()
            }
            rows = connection.execute(
                "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
            ).fetchall()
            recorded = {str(row[0]): (str(row[1]), int(row[2])) for row in rows}
            if set(recorded) != expected_tables or physical_tables != expected_tables:
                raise ValueError("InterPro table inventory does not match capabilities")
            for table_name, (role, row_count) in recorded.items():
                expected_role, expected_columns = _TABLE_CONTRACTS[table_name]
                if role != expected_role:
                    raise ValueError("InterPro table role inventory is unsupported")
                columns = [
                    (str(row[1]), str(row[2]))
                    for row in connection.execute(
                        f'PRAGMA table_info("{table_name}")'
                    ).fetchall()
                ]
                if columns != expected_columns:
                    raise ValueError(
                        f"InterPro table schema is unsupported: {table_name}"
                    )
                actual_count = connection.execute(
                    f'SELECT count(*) FROM "{table_name}"'
                ).fetchone()
                if actual_count is None or int(actual_count[0]) != row_count:
                    raise ValueError(f"InterPro row-count drift: {table_name}")

            observed_mappings = {
                tuple(str(value) for value in row)
                for row in connection.execute(
                    "SELECT table_name, source_column, output_column, reason "
                    "FROM _bioextract.column_mapping"
                ).fetchall()
            }
            expected_mappings = {
                (table, source, output, "generated_snake_case")
                for table in expected_tables
                for source, output in _source_output_columns(table)
                if source != output
            }
            if observed_mappings != expected_mappings:
                raise ValueError("InterPro column provenance inventory is unsupported")
            return capabilities, metadata.get("bioextract.release_version")
    except duckdb.Error as error:
        raise ValueError(f"Cannot open InterPro DuckDB publication: {path}") from error


def _source_output_columns(table: str) -> list[tuple[str, str]]:
    sources = {
        "mapping": list(_MAPPING_COLUMNS),
        "protein_term": ["UniProtId", "PfamId"],
        "term": ["PfamId", "PfamName"],
        "term_xref": ["PfamId", "InterProId", "InterProName", "InterProType"],
    }[table]
    outputs = [name for name, _type in _TABLE_CONTRACTS[table][1]]
    return list(zip(sources, outputs, strict=True))
