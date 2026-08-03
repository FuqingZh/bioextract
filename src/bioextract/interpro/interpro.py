from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._publication import DuckDBWriteResult, ParquetWriteResult
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
    InterProTidyConfig,
)
from .util import (
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
    file_protein2ipr: Path
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
    ) -> TidyDataset:
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
            enrichment fields null. The Pfam configuration requires both
            logical source roles and does not depend on a previously written
            canonical mapping. Paths never supply release identity; official
            XML metadata does.

        Examples:
            The XML below includes an INTERPRO release record whose official
            version is published as the release identity.

            Build the default mapping plan:

            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
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

        release_version = (
            read_interpro_release_version(self.snapshot.file_interpro_xml)
            if self.snapshot.file_interpro_xml is not None
            else None
        )
        if self.snapshot.file_interpro_xml is not None:
            validate_mapping_xml_relationships(
                self.snapshot.file_protein2ipr,
                df_interpro_entry=self.xml_frame("entry"),
                df_interpro_member=self.xml_frame("member"),
            )
        return TidyDataset(
            frames={
                "mapping": scan_mapping_frame(
                    self.snapshot.file_protein2ipr,
                    df_interpro_entry=self.xml_frame("entry"),
                    df_interpro_member=self.xml_frame("member"),
                )
            },
            source=self._tidy_sources(),
            resource_schema_version=SCHEMA_VERSION,
            source_schema_profile="interpro-protein2ipr-v1",
            build_id_prefix=f"interpro-mapping-{self.snapshot.file_protein2ipr.stem}",
            assets=tuple(
                TidyAsset(path=path, kind=kind, frame_name=frame_name)
                for path, kind, frame_name in ASSET_SPECS
            ),
            resource_name="interpro",
            release_version=release_version,
            release_version_source=(
                "official_metadata" if release_version is not None else None
            ),
        )

    def write_parquet(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: str = "fail",
    ) -> ParquetWriteResult:
        """Atomically publish the canonical InterPro mapping as one Parquet.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz"
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_parquet(
            ...         Path(dir_out) / "interpro.parquet"
            ...     )
            ...     result.resource_name.startswith("interpro-")
            True
        """
        return self.build_tidy().write_parquet(path, if_exists=if_exists)

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        config: InterProTidyConfig = "pfam",
        if_exists: str = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish a multi-relation InterPro product as DuckDB.

        Examples:
            >>> from tempfile import TemporaryDirectory
            >>> db = InterProDatabase.from_mapping_files(
            ...     protein_to_interpro="fixtures/interpro/108.0/raw/protein2ipr.dat.gz",
            ...     interpro_xml="fixtures/interpro/108.0/raw/interpro.xml.gz",
            ... )
            >>> with TemporaryDirectory() as dir_out:
            ...     result = db.write_duckdb(
            ...         Path(dir_out) / "interpro_pfam.duckdb"
            ...     )
            ...     result.tables
            ('protein_term', 'term', 'term_xref')
        """
        if config != "pfam":
            raise ValueError("write_duckdb() currently supports config='pfam' only")
        return self.build_tidy(config=config).write_duckdb(
            path,
            if_exists=if_exists,
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
                path=self.snapshot.file_protein2ipr,
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
            self._df_mapping = select_mapping_frame(
                self.dataset.snapshot.file_protein2ipr,
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
