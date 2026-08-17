from __future__ import annotations

import copy
import json
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import TYPE_CHECKING, Any, cast

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import (
    _RelationScanRequest,  # pyright: ignore[reportPrivateUsage]  # typed source boundary
    register_replayable_source,
)
from bioextract._publication import (
    METADATA_SCHEMA_VERSION,
    validate_duckdb_metadata_v2,
)

from ._knowledgebase import (
    RESOURCE_SCHEMA_VERSION,
    SOURCE_SCHEMA_PROFILE,
    TABLE_SCHEMAS,
)
from .constant import (
    COLS_IDMAPPING_SELECTED,
    IDMAPPING_SOURCE_SCHEMA_PROFILE,
    SCHEMA_MAPPING,
    SCHEMA_VERSION,
)
from .util import normalize_taxids

if TYPE_CHECKING:
    from .uniprot import UniProtDatabase

_NAMESPACES = {
    "uniprot",
    "entry_name",
    "gene_name",
    "gene_id",
    "refseq",
    "ensembl",
    "isoform_id",
}

# DuckDB can evaluate a bounded parameterized predicate more efficiently than
# building and hashing a temporary relation for the small selected-ID cases
# that motivated this API. Larger selections keep the temporary relation path
# so the physical strategy does not depend on an unbounded parameter list.
_SELECTED_MAPPING_IN_THRESHOLD = 1_000


@dataclass(frozen=True, slots=True)
class _ExtractSpec:
    relation_columns: str
    join: str
    order_by: str
    condition: str = ""
    parameters: tuple[object, ...] = ()


@dataclass(slots=True)
class _AnchorState:
    lock: Lock = field(default_factory=Lock)
    frame: pl.DataFrame | None = None


@dataclass(slots=True)
class _MappingAnchorState:
    lock: Lock = field(default_factory=Lock)
    matched_ids: frozenset[str] | None = None


@dataclass(frozen=True, slots=True)
class UniProtMappingSelection:
    """Selected-ID access over an idmapping source or publication."""

    database: UniProtDatabase = field(repr=False)
    input_ids: tuple[str, ...]
    taxon_ids: tuple[str, ...] = ()
    _anchor_state: _MappingAnchorState = field(
        default_factory=_MappingAnchorState,
        repr=False,
        compare=False,
    )

    def mappings(self) -> pl.LazyFrame:
        """Return only mapping rows whose UniProt ID was selected."""
        columns = ["input_id", "input_namespace", *COLS_IDMAPPING_SELECTED]
        schema: SchemaDict = {
            "input_id": pl.String,
            "input_namespace": pl.String,
            **SCHEMA_MAPPING,
        }
        if self.database._mode != "idmapping_publication":  # pyright: ignore[reportPrivateUsage]  # paired database boundary
            frame = self.database.scan_mapping(taxon_ids=self.taxon_ids)
            return (
                frame.filter(pl.col("uniprot_id").is_in(self.input_ids))
                .with_columns(
                    pl.col("uniprot_id").alias("input_id"),
                    pl.lit("uniprot").alias("input_namespace"),
                )
                .select(columns)
            )
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_selected_mapping_batches(
                snapshot,
                schema=schema,
                request=request,
            ),
            is_pure=True,
        )

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return selected UniProt IDs absent from the mapping scope."""
        schema: SchemaDict = {
            "input_id": pl.String,
            "input_namespace": pl.String,
            "reason": pl.String,
        }
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_selected_mapping_unmatched_batches(
                snapshot,
                schema=schema,
                request=request,
            ),
            is_pure=True,
        )


@dataclass(frozen=True, slots=True)
class UniProtSelection:
    """A reusable identifier selection with stable lazy domain relations.

    Examples:
        Resolve one accession, then extract its ordered accession relation:

        >>> selection = database.select_ids(  # doctest: +SKIP
        ...     ["P04637"], namespace="uniprot", taxon_ids=["9606"]
        ... )
        >>> selection.accessions().select(  # doctest: +SKIP
        ...     "primary_accession", "accession", "is_primary"
        ... )
        shape: (...)
    """

    database: UniProtDatabase = field(repr=False)
    namespace: str
    input_ids: tuple[str, ...]
    group_membership: tuple[tuple[str | None, str], ...]
    group_ids: tuple[str, ...] = ()
    taxon_ids: tuple[str, ...] = ()
    _is_grouped: bool = field(default=False, repr=False)
    _anchor_state: _AnchorState = field(
        default_factory=_AnchorState,
        repr=False,
        compare=False,
    )

    def proteins(self) -> pl.LazyFrame:
        """Return matched protein facts lazily.

        Examples:
            >>> selection.proteins().select(  # doctest: +SKIP
            ...     "primary_accession", "entry_name", "taxon_id"
            ... ).collect_schema().names()
            ['primary_accession', 'entry_name', 'taxon_id']
        """
        return self._lazy_relation(
            [
                "entry_name",
                "is_reviewed",
                "taxon_id",
                "protein_existence",
                "sequence_length",
                "molecular_weight",
                "sequence_version",
                "entry_version",
            ],
            _ExtractSpec(
                "p.entry_name AS entry_name, p.is_reviewed AS is_reviewed, "
                "p.taxon_id AS taxon_id, p.protein_existence AS protein_existence, "
                "p.sequence_length AS sequence_length, "
                "p.molecular_weight AS molecular_weight, "
                "p.sequence_version AS sequence_version, p.entry_version AS entry_version",
                "JOIN protein p USING (primary_accession)",
                "p.primary_accession",
            ),
        )

    def accessions(self) -> pl.LazyFrame:
        """Return ordered primary and secondary accessions lazily.

        Examples:
            >>> selection.accessions().select(  # doctest: +SKIP
            ...     "accession", "is_primary"
            ... ).collect_schema().names()
            ['accession', 'is_primary']
        """
        return self._lazy_relation(
            ["accession", "accession_order", "is_primary"],
            _ExtractSpec(
                "r.accession AS accession, r.accession_order AS accession_order, "
                "r.is_primary AS is_primary",
                "JOIN protein_accession r USING (primary_accession)",
                "r.accession_order, r.accession",
            ),
        )

    def protein_names(self) -> pl.LazyFrame:
        """Return official protein names lazily.

        Examples:
            >>> selection.protein_names().select(  # doctest: +SKIP
            ...     "name_type", "name"
            ... ).collect_schema().names()
            ['name_type', 'name']
        """
        return self._lazy_relation(
            ["name_type", "name", "name_order"],
            _ExtractSpec(
                "r.name_type AS name_type, r.name AS name, r.name_order AS name_order",
                "JOIN protein_name r USING (primary_accession)",
                "r.name_order, r.name",
            ),
        )

    def gene_names(self) -> pl.LazyFrame:
        """Return parsed gene names lazily.

        Examples:
            >>> selection.gene_names().select(  # doctest: +SKIP
            ...     "name_type", "name"
            ... ).collect_schema().names()
            ['name_type', 'name']
        """
        return self._lazy_relation(
            ["name_type", "name", "name_order"],
            _ExtractSpec(
                "r.name_type AS name_type, r.name AS name, r.name_order AS name_order",
                "JOIN gene_name r USING (primary_accession)",
                "r.name_order, r.name",
            ),
        )

    def ec_numbers(self) -> pl.LazyFrame:
        """Return distinct EC annotations lazily.

        Examples:
            >>> selection.ec_numbers().select("ec_number").collect_schema().names()  # doctest: +SKIP
            ['ec_number']
        """
        return self._lazy_relation(
            ["ec_number"],
            _ExtractSpec(
                "r.ec_number AS ec_number",
                "JOIN protein_ec_number r USING (primary_accession)",
                "r.ec_number",
            ),
        )

    def go_annotations(self) -> pl.LazyFrame:
        """Return GO annotations lazily.

        Examples:
            >>> selection.go_annotations().select(  # doctest: +SKIP
            ...     "go_id", "aspect", "evidence_code"
            ... ).collect_schema().names()
            ['go_id', 'aspect', 'evidence_code']
        """
        return self._lazy_relation(
            ["go_id", "aspect", "term_name", "evidence_code", "evidence_source"],
            _ExtractSpec(
                "r.go_id AS go_id, r.aspect AS aspect, r.term_name AS term_name, "
                "r.evidence_code AS evidence_code, r.evidence_source AS evidence_source",
                "JOIN protein_go_annotation r USING (primary_accession)",
                "r.go_id, r.aspect",
            ),
        )

    def cross_references(self, databases: Iterable[str] | None = None) -> pl.LazyFrame:
        """Return cross-references lazily, optionally filtered by database.

        Examples:
            >>> selection.cross_references(["PDB"]).select(  # doctest: +SKIP
            ...     "database", "external_id"
            ... ).collect_schema().names()
            ['database', 'external_id']
        """
        values = tuple(dict.fromkeys(databases or ()))
        condition = ""
        parameters: tuple[object, ...] = ()
        if values:
            condition = f" WHERE r.database IN ({','.join('?' for _ in values)})"
            parameters = values
        return self._lazy_relation(
            ["database", "external_id", "properties", "isoform_id"],
            _ExtractSpec(
                "r.database AS database, r.external_id AS external_id, "
                "r.properties AS properties, r.isoform_id AS isoform_id",
                "JOIN protein_cross_reference r USING (primary_accession)",
                "r.database, r.external_id, r.isoform_id",
                condition,
                parameters,
            ),
        )

    def comments(self, comment_types: Iterable[str] | None = None) -> pl.LazyFrame:
        """Return comments lazily, optionally filtered by comment type.

        Examples:
            >>> selection.comments(["FUNCTION"]).select(  # doctest: +SKIP
            ...     "comment_type", "comment_text"
            ... ).collect_schema().names()
            ['comment_type', 'comment_text']
        """
        values = tuple(dict.fromkeys(comment_types or ()))
        condition = ""
        parameters: tuple[object, ...] = ()
        if values:
            condition = f" WHERE r.comment_type IN ({','.join('?' for _ in values)})"
            parameters = values
        return self._lazy_relation(
            ["comment_id", "comment_type", "comment_text"],
            _ExtractSpec(
                "r.comment_id AS comment_id, r.comment_type AS comment_type, "
                "r.comment_text AS comment_text",
                "JOIN protein_comment r USING (primary_accession)",
                "r.comment_id",
                condition,
                parameters,
            ),
        )

    def subcellular_locations(self) -> pl.LazyFrame:
        """Return individually parsed locations lazily.

        Examples:
            >>> selection.subcellular_locations().select(  # doctest: +SKIP
            ...     "location", "note"
            ... ).collect_schema().names()
            ['location', 'note']
        """
        return self._lazy_relation(
            ["location", "note"],
            _ExtractSpec(
                "r.location AS location, r.note AS note",
                "JOIN protein_subcellular_location r USING (primary_accession)",
                "r.comment_id, r.location",
            ),
        )

    def keywords(self) -> pl.LazyFrame:
        """Return keywords lazily in source-defined order.

        Examples:
            >>> selection.keywords().select(  # doctest: +SKIP
            ...     "keyword", "keyword_order"
            ... ).collect_schema().names()
            ['keyword', 'keyword_order']
        """
        return self._lazy_relation(
            ["keyword", "keyword_order"],
            _ExtractSpec(
                "r.keyword AS keyword, r.keyword_order AS keyword_order",
                "JOIN protein_keyword r USING (primary_accession)",
                "r.keyword_order, r.keyword",
            ),
        )

    def sequences(self, sequence_type: str = "canonical") -> pl.LazyFrame:
        """Return canonical, isoform, or all sequences lazily.

        Examples:
            >>> selection.sequences("all").select(  # doctest: +SKIP
            ...     "sequence_id", "sequence_type", "sha256"
            ... ).collect_schema().names()
            ['sequence_id', 'sequence_type', 'sha256']
        """
        if sequence_type not in {"canonical", "isoform", "all"}:
            raise ValueError("sequence_type must be canonical, isoform, or all")
        condition = "" if sequence_type == "all" else " WHERE r.sequence_type = ?"
        parameters: tuple[object, ...] = (
            () if sequence_type == "all" else (sequence_type,)
        )
        return self._lazy_relation(
            ["sequence_id", "sequence_type", "sequence", "length", "crc64", "sha256"],
            _ExtractSpec(
                "r.sequence_id AS sequence_id, r.sequence_type AS sequence_type, "
                "r.sequence AS sequence, r.length AS length, r.crc64 AS crc64, "
                "r.sha256 AS sha256",
                "JOIN protein_sequence r USING (primary_accession)",
                "CASE r.sequence_type WHEN 'canonical' THEN 0 ELSE 1 END, r.sequence_id",
                condition,
                parameters,
            ),
        )

    def isoforms(self) -> pl.LazyFrame:
        """Return isoform definitions lazily.

        Examples:
            >>> selection.isoforms().select(  # doctest: +SKIP
            ...     "isoform_id", "sequence_status", "sequence_id"
            ... ).collect_schema().names()
            ['isoform_id', 'sequence_status', 'sequence_id']
        """
        return self._lazy_relation(
            ["isoform_id", "name", "isoform_order", "sequence_status", "sequence_id"],
            _ExtractSpec(
                "r.isoform_id AS isoform_id, r.name AS name, "
                "r.isoform_order AS isoform_order, "
                "r.sequence_status AS sequence_status, r.sequence_id AS sequence_id",
                "JOIN protein_isoform r USING (primary_accession)",
                "r.isoform_order, r.isoform_id",
            ),
        )

    def isoform_identifiers(self) -> pl.LazyFrame:
        """Return every current or old official IsoId lazily.

        Examples:
            >>> selection.isoform_identifiers().select(  # doctest: +SKIP
            ...     "isoform_id", "identifier", "is_main"
            ... ).collect_schema().names()
            ['isoform_id', 'identifier', 'is_main']
        """
        return self._lazy_relation(
            ["isoform_id", "identifier", "identifier_order", "is_main"],
            _ExtractSpec(
                "r.isoform_id AS isoform_id, r.identifier AS identifier, "
                "r.identifier_order AS identifier_order, r.is_main AS is_main",
                "JOIN protein_isoform_identifier r USING (primary_accession)",
                "r.isoform_id, r.identifier_order",
            ),
        )

    def sequence_variations(self) -> pl.LazyFrame:
        """Return DAT VAR_SEQ features lazily.

        Examples:
            >>> selection.sequence_variations().select(  # doctest: +SKIP
            ...     "variation_id", "start_position", "end_position"
            ... ).collect_schema().names()
            ['variation_id', 'start_position', 'end_position']
        """
        return self._lazy_relation(
            ["variation_id", "start_position", "end_position", "note"],
            _ExtractSpec(
                "r.variation_id AS variation_id, r.start_position AS start_position, "
                "r.end_position AS end_position, r.note AS note",
                "JOIN protein_sequence_variation r USING (primary_accession)",
                "r.start_position, r.end_position, r.variation_id",
            ),
        )

    def isoform_variations(self) -> pl.LazyFrame:
        """Return ordered isoform-to-VSP relationships lazily.

        Examples:
            >>> selection.isoform_variations().select(  # doctest: +SKIP
            ...     "isoform_id", "variation_id", "variation_order"
            ... ).collect_schema().names()
            ['isoform_id', 'variation_id', 'variation_order']
        """
        return self._lazy_relation(
            ["isoform_id", "variation_id", "variation_order"],
            _ExtractSpec(
                "r.isoform_id AS isoform_id, r.variation_id AS variation_id, "
                "r.variation_order AS variation_order",
                "JOIN protein_isoform_variation r USING (primary_accession)",
                "r.isoform_id, r.variation_order, r.variation_id",
            ),
        )

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return requested identifiers that matched no protein lazily.

        Examples:
            >>> selection.unmatched_ids().select(  # doctest: +SKIP
            ...     "input_id", "input_namespace", "reason"
            ... ).collect_schema().names()
            ['input_id', 'input_namespace', 'reason']
        """
        columns = ["input_id", "input_namespace", "reason"]
        if self._is_grouped:
            columns.insert(0, "group_id")
        snapshot = copy.copy(self)
        return register_replayable_source(
            schema=dict.fromkeys(columns, pl.String),
            batches=lambda request: _iter_unmatched_batches(
                snapshot,
                request=request,
            ),
            is_pure=True,
        )

    def _lazy_relation(
        self,
        relation_columns: Iterable[str],
        spec: _ExtractSpec,
    ) -> pl.LazyFrame:
        snapshot = copy.copy(self)
        columns = ["input_id", "input_namespace", "primary_accession"]
        if self._is_grouped:
            columns.insert(0, "group_id")
        schema = self._schema_for_relation([*columns, *relation_columns])
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_relation_batches(
                snapshot,
                spec=spec,
                schema=schema,
                request=request,
            ),
            is_pure=True,
        )

    @staticmethod
    def _schema_for_relation(columns: Iterable[str]) -> SchemaDict:
        """Return stable Polars types for selection and domain columns."""
        relation_types: dict[str, Any] = {
            "group_id": pl.String,
            "input_id": pl.String,
            "input_namespace": pl.String,
            "primary_accession": pl.String,
            "reason": pl.String,
        }
        for table_schema in TABLE_SCHEMAS.values():
            for name, dtype in table_schema.items():
                if name not in relation_types:
                    relation_types[name] = dtype
        return {column: relation_types[column] for column in columns}

    def _matches(self) -> pl.DataFrame:
        matches = self._identifier_matches()
        membership = pl.DataFrame(
            self.group_membership,
            schema={"group_id": pl.String, "input_id": pl.String},
            orient="row",
        )
        if membership.is_empty() or matches.is_empty():
            return pl.DataFrame(
                schema={
                    "group_id": pl.String,
                    "input_id": pl.String,
                    "input_namespace": pl.String,
                    "primary_accession": pl.String,
                }
            )
        return (
            membership.join(matches, on="input_id", how="inner")
            .select("group_id", "input_id", "input_namespace", "primary_accession")
            .sort("group_id", "input_id", "primary_accession")
        )

    def _identifier_matches(self) -> pl.DataFrame:
        with self._anchor_state.lock:
            if self._anchor_state.frame is None:
                self._anchor_state.frame = self._query_identifier_matches()
            return self._anchor_state.frame.clone()

    def _query_identifier_matches(self) -> pl.DataFrame:
        inputs = pl.DataFrame(
            {"input_id": self.input_ids},
            schema={"input_id": pl.String},
        )
        with self.database.connect() as connection:
            connection.execute("CREATE TEMP TABLE _inputs(input_id VARCHAR)")
            if inputs.height:
                connection.executemany("INSERT INTO _inputs VALUES (?)", inputs.rows())
            taxon_filter = ""
            parameters: list[object] = [self.namespace, self.namespace]
            if self.taxon_ids:
                taxon_filter = (
                    f" WHERE pr.taxon_id IN ({','.join('?' for _ in self.taxon_ids)})"
                )
                parameters.extend(self.taxon_ids)
            cursor = connection.execute(
                f"""
                SELECT DISTINCT i.input_id, ? AS input_namespace,
                       p.primary_accession
                FROM _inputs i
                JOIN protein_identifier p
                  ON p.identifier = i.input_id AND p.namespace = ?
                JOIN protein pr USING (primary_accession)
                {taxon_filter}
                ORDER BY i.input_id, p.primary_accession
                """,
                parameters,
            )
            return _cursor_frame(cursor)


def validate_publication(path: Path) -> tuple[str, Mapping[str, str]]:
    with duckdb.connect(str(path), read_only=True) as connection:
        metadata_tables = {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema='_bioextract'"
            ).fetchall()
        }
        required_metadata_tables = {
            "metadata",
            "source_file",
            "table_info",
            "column_mapping",
            "validation_issue",
        }
        missing_metadata_tables = sorted(required_metadata_tables - metadata_tables)
        if missing_metadata_tables:
            raise ValueError(
                "UniProtKB publication lacks metadata tables: "
                f"{missing_metadata_tables}"
            )
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        if (
            metadata.get("bioextract.metadata_schema_version")
            != METADATA_SCHEMA_VERSION
        ):
            raise ValueError("Unsupported UniProt metadata schema version")
        required_v2 = {
            "bioextract.resource_name",
            "bioextract.resource_schema_version",
            "bioextract.source_schema_profile",
            "bioextract.package_version",
            "bioextract.generated_at",
            "bioextract.validation_status",
            "bioextract.validation_issue_count",
        }
        missing_v2 = sorted(required_v2 - set(metadata))
        if missing_v2:
            raise ValueError(f"UniProt metadata v2 is missing keys: {missing_v2}")
        validate_duckdb_metadata_v2(connection, metadata)
        issue_count = connection.execute(
            "SELECT count(*) FROM _bioextract.validation_issue"
        ).fetchone()
        if issue_count is None or int(
            metadata.get("bioextract.validation_issue_count", "-1")
        ) != int(issue_count[0]):
            raise ValueError("UniProt validation issue count mismatch")
        if metadata.get("bioextract.resource_name") != "uniprot":
            raise ValueError("DuckDB file is not a bioextract UniProt publication")
        profile_metadata = (
            metadata.get("bioextract.resource_schema_version"),
            metadata.get("bioextract.source_schema_profile"),
        )
        if profile_metadata == (SCHEMA_VERSION, IDMAPPING_SOURCE_SCHEMA_PROFILE):
            profile = "idmapping"
            expected_schemas = {"mapping": SCHEMA_MAPPING}
            required_capabilities = {"bioextract.capability.mapping": "true"}
            expected_capability_keys = set(required_capabilities)
        elif profile_metadata == (RESOURCE_SCHEMA_VERSION, SOURCE_SCHEMA_PROFILE):
            profile = "knowledgebase"
            expected_schemas = TABLE_SCHEMAS
            required_capabilities = {
                "bioextract.capability.canonical_sequences": "true",
                "bioextract.capability.isoform_definitions": "true",
            }
            expected_capability_keys = {
                *required_capabilities,
                "bioextract.capability.isoform_sequences",
            }
        else:
            raise ValueError(
                "Unsupported UniProt source schema profile or resource schema"
            )
        capability_metadata = {
            key: value
            for key, value in metadata.items()
            if key.startswith("bioextract.capability.")
        }
        if any(
            capability_metadata.get(key) != value
            for key, value in required_capabilities.items()
        ):
            raise ValueError("Unsupported UniProt capability metadata")
        if profile == "idmapping" and capability_metadata != required_capabilities:
            raise ValueError("Unsupported UniProt idmapping capability inventory")
        if profile == "idmapping":
            source_roles = {
                str(row[0])
                for row in connection.execute(
                    "SELECT logical_name FROM _bioextract.source_file"
                ).fetchall()
            }
            if source_roles != {"idmapping_selected"}:
                raise ValueError("Unsupported UniProt idmapping source inventory")
            column_mapping_count = connection.execute(
                "SELECT count(*) FROM _bioextract.column_mapping"
            ).fetchone()
            if column_mapping_count != (0,):
                raise ValueError("Unsupported UniProt idmapping column provenance")
            try:
                scope = cast(object, json.loads(metadata["bioextract.scope"]))
            except (KeyError, json.JSONDecodeError) as error:
                raise ValueError("Unsupported UniProt idmapping taxon scope") from error
            scope_mapping = (
                cast(dict[str, object], scope) if isinstance(scope, dict) else {}
            )
            valid_all_taxa = (
                set(scope_mapping) == {"all_taxa"}
                and isinstance(scope_mapping["all_taxa"], bool)
                and scope_mapping["all_taxa"] is True
            )
            taxon_ids_value = scope_mapping.get("taxon_ids")
            taxon_ids = (
                cast(list[object], taxon_ids_value)
                if isinstance(taxon_ids_value, list)
                else None
            )
            valid_taxon_ids = (
                taxon_ids is not None
                and bool(taxon_ids)
                and all(
                    isinstance(taxon_id, str)
                    and bool(taxon_id)
                    and taxon_id == taxon_id.strip()
                    for taxon_id in taxon_ids
                )
                and len(taxon_ids) == len(set(taxon_ids))
                and set(scope_mapping) == {"taxon_ids"}
            )
            if not valid_all_taxa and not valid_taxon_ids:
                raise ValueError("Unsupported UniProt idmapping taxon scope")
        if profile == "knowledgebase" and (
            set(capability_metadata) != expected_capability_keys
            or capability_metadata["bioextract.capability.isoform_sequences"]
            not in {"true", "false"}
        ):
            raise ValueError("Unsupported UniProtKB capability inventory")
        relations = {
            (str(row[0]), str(row[1]))
            for row in connection.execute(
                "SELECT table_name, table_type FROM information_schema.tables "
                "WHERE table_schema='main'"
            ).fetchall()
        }
        expected_tables = set(expected_schemas)
        if relations != {(table, "BASE TABLE") for table in expected_tables}:
            raise ValueError(
                "UniProtKB relation inventory mismatch: "
                f"expected={sorted(expected_tables)}, actual={sorted(relations)}"
            )
        recorded_rows = connection.execute(
            "SELECT table_name, table_role, row_count FROM _bioextract.table_info"
        ).fetchall()
        recorded_tables = [str(row[0]) for row in recorded_rows]
        if len(recorded_tables) != len(set(recorded_tables)):
            raise ValueError("UniProtKB table_info contains duplicate relations")
        if set(recorded_tables) != expected_tables:
            raise ValueError(
                "UniProtKB table_info inventory mismatch: "
                f"expected={sorted(expected_tables)}, actual={sorted(recorded_tables)}"
            )
        if any(str(row[1]) != "canonical" for row in recorded_rows):
            raise ValueError("UniProt table_info role inventory is unsupported")
        for table, schema in expected_schemas.items():
            actual_columns = connection.execute(
                """
                SELECT column_name, data_type, is_nullable
                FROM information_schema.columns
                WHERE table_schema = 'main' AND table_name = ?
                ORDER BY ordinal_position
                """,
                [table],
            ).fetchall()
            expected_columns = [
                (column, _duckdb_type(dtype), "YES") for column, dtype in schema.items()
            ]
            if actual_columns != expected_columns:
                raise ValueError(
                    f"UniProtKB physical schema mismatch for {table}: "
                    f"expected={expected_columns}, actual={actual_columns}"
                )
        return profile, metadata


def make_selection(
    database: UniProtDatabase,
    ids: Iterable[str],
    *,
    namespace: str,
    groups: Mapping[str, Iterable[str]] | None = None,
    taxon_ids: Iterable[str | int] | None = None,
) -> UniProtSelection:
    if namespace not in _NAMESPACES:
        raise ValueError(f"Unsupported UniProt identifier namespace: {namespace}")
    if groups is None:
        input_ids = _normalize_ids(ids)
        membership = tuple((None, value) for value in input_ids)
        group_ids: tuple[str, ...] = ()
    else:
        normalized_groups: dict[str, Iterable[str]] = {}
        for group, values in groups.items():
            normalized_group = str(group).strip()
            if not normalized_group:
                raise ValueError("UniProt group labels must be non-empty")
            if normalized_group in normalized_groups:
                raise ValueError(
                    "UniProt group labels must be unique after normalization"
                )
            normalized_groups[normalized_group] = values
        group_ids = tuple(sorted(normalized_groups))
        membership = tuple(
            sorted(
                (group, value)
                for group, values in normalized_groups.items()
                for value in _normalize_ids(values)
            )
        )
        input_ids = tuple(sorted({value for _, value in membership}))
    normalized_taxa = normalize_taxids(tuple(taxon_ids or ()))
    return UniProtSelection(
        database=database,
        namespace=namespace,
        input_ids=input_ids,
        group_membership=membership,
        group_ids=group_ids,
        taxon_ids=normalized_taxa,
        _is_grouped=groups is not None,
    )


def _normalize_ids(ids: Iterable[str]) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            str(identifier).strip() for identifier in ids if str(identifier).strip()
        )
    )


def _duckdb_type(dtype: object) -> str:
    if dtype == pl.String:
        return "VARCHAR"
    if dtype == pl.Int64:
        return "BIGINT"
    if dtype == pl.Boolean:
        return "BOOLEAN"
    raise TypeError(f"Unsupported UniProtKB physical type: {dtype}")


def _cursor_frame(cursor: duckdb.DuckDBPyConnection) -> pl.DataFrame:
    schema = {
        description[0]: _polars_type(str(description[1]))
        for description in cursor.description
    }
    rows = cursor.fetchall()
    return pl.DataFrame(rows, schema=schema, orient="row")


def _iter_selected_mapping_batches(
    selection: UniProtMappingSelection,
    *,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    if not selection.input_ids:
        return
    connection = selection.database.connect()  # pyright: ignore[reportPrivateUsage]  # selected mapping boundary
    reader: Any = None
    try:
        requested = _requested_columns(request.columns, schema)
        output_columns = list(schema) if requested is None else requested
        source_sql, parameters, input_expression = _selected_mapping_source(
            connection,
            selection,
        )
        projection: list[str] = []
        for name in output_columns:
            if name == "input_id":
                projection.append(f'{input_expression} AS "input_id"')
            elif name == "input_namespace":
                projection.append("'uniprot' AS input_namespace")
            else:
                projection.append(f'mapping."{name}" AS "{name}"')
        query = f"""
            SELECT {", ".join(projection)}
            {source_sql}
        """
        query += f" ORDER BY {input_expression}"
        result = connection.execute(query, parameters)
        reader = _uniprot_arrow_reader(result, request.effective_batch_size)
        for record_batch in reader:
            frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
            yield frame.cast(  # type: ignore[reportArgumentType]
                {name: schema[name] for name in frame.columns}, strict=False
            )
    finally:
        if reader is not None:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
        connection.close()


def _iter_selected_mapping_unmatched_batches(
    selection: UniProtMappingSelection,
    *,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    if not selection.input_ids:
        return
    if selection.database._mode == "idmapping_publication":  # pyright: ignore[reportPrivateUsage]  # paired database boundary
        matched = set(_publication_selected_mapping_ids(selection))
    else:
        matched = set(
            selection.database.scan_mapping(taxon_ids=selection.taxon_ids)
            .filter(pl.col("uniprot_id").is_in(selection.input_ids))
            .select("uniprot_id")
            .unique()
            .collect()
            .get_column("uniprot_id")
            .to_list()
        )
    frame = pl.DataFrame(
        {
            "input_id": [
                value for value in selection.input_ids if value not in matched
            ],
            "input_namespace": [
                "uniprot" for value in selection.input_ids if value not in matched
            ],
            "reason": [
                "not_found" for value in selection.input_ids if value not in matched
            ],
        },
        schema=schema,
    )
    requested = _requested_columns(request.columns, schema)
    if requested is not None:
        frame = frame.select(requested)
    yield from _frame_slices(frame, request.effective_batch_size)


def _selected_mapping_source(
    connection: duckdb.DuckDBPyConnection,
    selection: UniProtMappingSelection,
) -> tuple[str, list[str], str]:
    parameters: list[str]
    if len(selection.input_ids) <= _SELECTED_MAPPING_IN_THRESHOLD:
        source_sql = (
            "FROM mapping AS mapping WHERE mapping.uniprot_id IN ("
            + ",".join("?" for _ in selection.input_ids)
            + ")"
        )
        parameters = list(selection.input_ids)
        input_expression = "mapping.uniprot_id"
    else:
        connection.execute(
            "CREATE TEMP TABLE _selected_uniprot(input_id VARCHAR PRIMARY KEY)"
        )
        connection.executemany(
            "INSERT INTO _selected_uniprot VALUES (?)",
            [(value,) for value in selection.input_ids],
        )
        source_sql = (
            "FROM mapping AS mapping JOIN _selected_uniprot AS selected "
            "ON mapping.uniprot_id = selected.input_id"
        )
        parameters = []
        input_expression = "selected.input_id"
    if selection.taxon_ids:
        source_sql += " AND" if " WHERE " in source_sql else " WHERE"
        source_sql += (
            " mapping.tax_id IN (" + ",".join("?" for _ in selection.taxon_ids) + ")"
        )
        parameters.extend(selection.taxon_ids)
    return source_sql, parameters, input_expression


def _publication_selected_mapping_ids(
    selection: UniProtMappingSelection,
) -> frozenset[str]:
    state = selection._anchor_state  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    with state.lock:
        if state.matched_ids is not None:
            return state.matched_ids
        connection = selection.database.connect()
        try:
            source_sql, parameters, input_expression = _selected_mapping_source(
                connection,
                selection,
            )
            matched = frozenset(
                str(row[0])
                for row in connection.execute(
                    f"SELECT DISTINCT {input_expression} {source_sql}",
                    parameters,
                ).fetchall()
            )
        finally:
            connection.close()
        state.matched_ids = matched
        return matched


def _iter_relation_batches(
    selection: UniProtSelection,
    *,
    spec: _ExtractSpec,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    matches = selection._matches()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    if matches.is_empty():
        return

    connection = selection.database.connect()
    reader: Any = None
    try:
        connection.execute(
            "CREATE TEMP TABLE _selection("
            "group_id VARCHAR, input_id VARCHAR, input_namespace VARCHAR, "
            "primary_accession VARCHAR)"
        )
        connection.executemany(
            "INSERT INTO _selection VALUES (?, ?, ?, ?)",
            matches.rows(),
        )
        inner_query = (
            "SELECT DISTINCT s.group_id, s.input_id, s.input_namespace, "
            "s.primary_accession AS primary_accession, "
            f"{spec.relation_columns} FROM _selection s {spec.join}{spec.condition} "
            "ORDER BY s.group_id, s.input_id, s.primary_accession, "
            f"{spec.order_by}"
        )
        requested = _requested_columns(request.columns, schema)
        query = (
            inner_query if requested is None else _project_query(inner_query, requested)
        )
        result = connection.execute(query, list(spec.parameters))
        reader = _uniprot_arrow_reader(result, request.effective_batch_size)
        for record_batch in reader:
            frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
            frame = frame.cast(  # type: ignore[reportArgumentType]
                {name: schema[name] for name in frame.columns if name in schema},
                strict=False,
            )
            if not selection._is_grouped and "group_id" in frame.columns:  # pyright: ignore[reportPrivateUsage]
                frame = frame.drop("group_id")
            yield frame
    finally:
        if reader is not None:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
        connection.close()


def _iter_unmatched_batches(
    selection: UniProtSelection,
    *,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    matched = set(
        selection._identifier_matches().get_column("input_id").to_list()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    )
    rows = [
        {
            "group_id": group,
            "input_id": input_id,
            "input_namespace": selection.namespace,
            "reason": "not_found",
        }
        for group, input_id in selection.group_membership
        if input_id not in matched
    ]
    schema: SchemaDict = {
        "group_id": pl.String,
        "input_id": pl.String,
        "input_namespace": pl.String,
        "reason": pl.String,
    }
    frame = pl.DataFrame(rows, schema=schema)
    if not selection._is_grouped:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        frame = frame.drop("group_id")
        schema = {name: dtype for name, dtype in schema.items() if name != "group_id"}
    requested = _requested_columns(request.columns, schema)
    if requested is not None:
        frame = frame.select(requested)
    yield from _frame_slices(frame, request.effective_batch_size)


def _requested_columns(
    columns: tuple[str, ...] | None,
    schema: SchemaDict,
) -> list[str] | None:
    if columns is None:
        return None
    selected = [name for name in columns if name in schema]
    return selected or None


def _project_query(query: str, columns: list[str]) -> str:
    projection = ", ".join(f'"{column}"' for column in columns)
    return f"SELECT {projection} FROM ({query}) AS _bioextract_relation"


def _frame_slices(frame: pl.DataFrame, batch_size: int) -> Iterator[pl.DataFrame]:
    for offset in range(0, frame.height, batch_size):
        yield frame.slice(offset, batch_size)


def _uniprot_arrow_reader(result: Any, batch_size: int) -> Any:
    to_arrow_reader = getattr(result, "to_arrow_reader", None)
    if to_arrow_reader is not None:
        return to_arrow_reader(batch_size)
    return result.fetch_record_batch(rows_per_batch=batch_size)


def _polars_type(type_name: str) -> type[pl.DataType]:
    if type_name == "BOOLEAN":
        return pl.Boolean
    if type_name in {"BIGINT", "INTEGER", "UBIGINT"}:
        return pl.Int64
    return pl.String
