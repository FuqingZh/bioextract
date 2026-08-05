from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, cast

import duckdb
import polars as pl

from bioextract._publication import validate_duckdb_metadata_v1

from ._knowledgebase import (
    RESOURCE_SCHEMA_VERSION,
    SOURCE_SCHEMA_PROFILE,
    TABLE_SCHEMAS,
)
from .constant import IDMAPPING_SOURCE_SCHEMA_PROFILE, SCHEMA_MAPPING, SCHEMA_VERSION
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


@dataclass(frozen=True, slots=True)
class UniProtSelection:
    """A reusable identifier selection with stable domain extractors.

    Examples:
        Resolve one accession, then extract its ordered accession relation:

        >>> selection = database.select_ids(  # doctest: +SKIP
        ...     ["P04637"], namespace="uniprot", taxon_ids=["9606"]
        ... )
        >>> selection.extract_accessions().select(  # doctest: +SKIP
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
    _matches_cache: pl.DataFrame | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def extract_proteins(self) -> pl.DataFrame:
        """Return matched protein facts in stable selection/accession order.

        Examples:
            >>> selection.extract_proteins().select(  # doctest: +SKIP
            ...     "primary_accession", "entry_name", "taxon_id"
            ... ).columns
            ['primary_accession', 'entry_name', 'taxon_id']
        """
        return self._extract(
            "protein",
            "p.entry_name AS entry_name, p.is_reviewed AS is_reviewed, "
            "p.taxon_id AS taxon_id, p.protein_existence AS protein_existence, "
            "p.sequence_length AS sequence_length, "
            "p.molecular_weight AS molecular_weight, "
            "p.sequence_version AS sequence_version, "
            "p.entry_version AS entry_version",
            "JOIN protein p USING (primary_accession)",
            order_by="p.primary_accession",
        )

    def extract_accessions(self) -> pl.DataFrame:
        """Return ordered primary and secondary accessions for every match.

        Examples:
            >>> selection.extract_accessions().select(  # doctest: +SKIP
            ...     "accession", "is_primary"
            ... ).columns
            ['accession', 'is_primary']
        """
        return self._extract(
            "protein_accession",
            "r.accession AS accession, r.accession_order AS accession_order, "
            "r.is_primary AS is_primary",
            "JOIN protein_accession r USING (primary_accession)",
            order_by="r.accession_order, r.accession",
        )

    def extract_protein_names(self) -> pl.DataFrame:
        """Return official protein names in source-defined name order.

        Examples:
            >>> selection.extract_protein_names().select(  # doctest: +SKIP
            ...     "name_type", "name"
            ... ).columns
            ['name_type', 'name']
        """
        return self._extract(
            "protein_name",
            "r.name_type AS name_type, r.name AS name, r.name_order AS name_order",
            "JOIN protein_name r USING (primary_accession)",
            order_by="r.name_order, r.name",
        )

    def extract_gene_names(self) -> pl.DataFrame:
        """Return parsed gene names in source-defined name order.

        Examples:
            >>> selection.extract_gene_names().select(  # doctest: +SKIP
            ...     "name_type", "name"
            ... ).columns
            ['name_type', 'name']
        """
        return self._extract(
            "gene_name",
            "r.name_type AS name_type, r.name AS name, r.name_order AS name_order",
            "JOIN gene_name r USING (primary_accession)",
            order_by="r.name_order, r.name",
        )

    def extract_ec_numbers(self) -> pl.DataFrame:
        """Return distinct EC annotations ordered by EC number.

        Examples:
            >>> selection.extract_ec_numbers().select("ec_number").columns  # doctest: +SKIP
            ['ec_number']
        """
        return self._extract(
            "protein_ec_number",
            "r.ec_number AS ec_number",
            "JOIN protein_ec_number r USING (primary_accession)",
            order_by="r.ec_number",
        )

    def extract_go_annotations(self) -> pl.DataFrame:
        """Return GO annotations ordered by GO identifier and aspect.

        Examples:
            >>> selection.extract_go_annotations().select(  # doctest: +SKIP
            ...     "go_id", "aspect", "evidence_code"
            ... ).columns
            ['go_id', 'aspect', 'evidence_code']
        """
        return self._extract(
            "protein_go_annotation",
            "r.go_id AS go_id, r.aspect AS aspect, r.term_name AS term_name, "
            "r.evidence_code AS evidence_code, r.evidence_source AS evidence_source",
            "JOIN protein_go_annotation r USING (primary_accession)",
            order_by="r.go_id, r.aspect",
        )

    def extract_cross_references(
        self, databases: Iterable[str] | None = None
    ) -> pl.DataFrame:
        """Return cross-references, optionally limited to database names.

        Examples:
            >>> selection.extract_cross_references(["PDB"]).select(  # doctest: +SKIP
            ...     "database", "external_id"
            ... ).columns
            ['database', 'external_id']
        """
        values = tuple(dict.fromkeys(databases or ()))
        condition = ""
        parameters: list[object] = []
        if values:
            condition = f" WHERE r.database IN ({','.join('?' for _ in values)})"
            parameters.extend(values)
        return self._extract(
            "protein_cross_reference",
            "r.database AS database, r.external_id AS external_id, "
            "r.properties AS properties, r.isoform_id AS isoform_id",
            "JOIN protein_cross_reference r USING (primary_accession)",
            condition=condition,
            parameters=parameters,
            order_by="r.database, r.external_id, r.isoform_id",
        )

    def extract_comments(
        self, comment_types: Iterable[str] | None = None
    ) -> pl.DataFrame:
        """Return comments, optionally limited to exact comment types.

        Examples:
            >>> selection.extract_comments(["FUNCTION"]).select(  # doctest: +SKIP
            ...     "comment_type", "comment_text"
            ... ).columns
            ['comment_type', 'comment_text']
        """
        values = tuple(dict.fromkeys(comment_types or ()))
        condition = ""
        parameters: list[object] = []
        if values:
            condition = f" WHERE r.comment_type IN ({','.join('?' for _ in values)})"
            parameters.extend(values)
        return self._extract(
            "protein_comment",
            "r.comment_id AS comment_id, r.comment_type AS comment_type, "
            "r.comment_text AS comment_text",
            "JOIN protein_comment r USING (primary_accession)",
            condition=condition,
            parameters=parameters,
            order_by="r.comment_id",
        )

    def extract_subcellular_locations(self) -> pl.DataFrame:
        """Return individually parsed locations and their optional notes.

        Examples:
            >>> selection.extract_subcellular_locations().select(  # doctest: +SKIP
            ...     "location", "note"
            ... ).columns
            ['location', 'note']
        """
        return self._extract(
            "protein_subcellular_location",
            "r.location AS location, r.note AS note",
            "JOIN protein_subcellular_location r USING (primary_accession)",
            order_by="r.comment_id, r.location",
        )

    def extract_keywords(self) -> pl.DataFrame:
        """Return keywords in source-defined keyword order.

        Examples:
            >>> selection.extract_keywords().select(  # doctest: +SKIP
            ...     "keyword", "keyword_order"
            ... ).columns
            ['keyword', 'keyword_order']
        """
        return self._extract(
            "protein_keyword",
            "r.keyword AS keyword, r.keyword_order AS keyword_order",
            "JOIN protein_keyword r USING (primary_accession)",
            order_by="r.keyword_order, r.keyword",
        )

    def extract_sequences(self, sequence_type: str = "canonical") -> pl.DataFrame:
        """Return canonical, isoform, or all sequences in stable type order.

        Examples:
            >>> selection.extract_sequences("all").select(  # doctest: +SKIP
            ...     "sequence_id", "sequence_type", "sha256"
            ... ).columns
            ['sequence_id', 'sequence_type', 'sha256']
        """
        if sequence_type not in {"canonical", "isoform", "all"}:
            raise ValueError("sequence_type must be canonical, isoform, or all")
        condition = "" if sequence_type == "all" else " WHERE r.sequence_type = ?"
        parameters: list[object] = [] if sequence_type == "all" else [sequence_type]
        return self._extract(
            "protein_sequence",
            "r.sequence_id AS sequence_id, r.sequence_type AS sequence_type, "
            "r.sequence AS sequence, r.length AS length, r.crc64 AS crc64, "
            "r.sha256 AS sha256",
            "JOIN protein_sequence r USING (primary_accession)",
            condition=condition,
            parameters=parameters,
            order_by=(
                "CASE r.sequence_type WHEN 'canonical' THEN 0 ELSE 1 END, r.sequence_id"
            ),
        )

    def extract_isoforms(self) -> pl.DataFrame:
        """Return isoform definitions with normalized status and sequence link.

        Examples:
            >>> selection.extract_isoforms().select(  # doctest: +SKIP
            ...     "isoform_id", "sequence_status", "sequence_id"
            ... ).columns
            ['isoform_id', 'sequence_status', 'sequence_id']
        """
        return self._extract(
            "protein_isoform",
            "r.isoform_id AS isoform_id, r.name AS name, "
            "r.isoform_order AS isoform_order, "
            "r.sequence_status AS sequence_status, r.sequence_id AS sequence_id",
            "JOIN protein_isoform r USING (primary_accession)",
            order_by="r.isoform_order, r.isoform_id",
        )

    def extract_isoform_identifiers(self) -> pl.DataFrame:
        """Return every current or old official IsoId for each product.

        Examples:
            >>> selection.extract_isoform_identifiers().select(  # doctest: +SKIP
            ...     "isoform_id", "identifier", "is_main"
            ... ).columns
            ['isoform_id', 'identifier', 'is_main']
        """
        return self._extract(
            "protein_isoform_identifier",
            "r.isoform_id AS isoform_id, r.identifier AS identifier, "
            "r.identifier_order AS identifier_order, r.is_main AS is_main",
            "JOIN protein_isoform_identifier r USING (primary_accession)",
            order_by="r.isoform_id, r.identifier_order",
        )

    def extract_sequence_variations(self) -> pl.DataFrame:
        """Return DAT VAR_SEQ features ordered by coordinates and VSP ID.

        Examples:
            >>> selection.extract_sequence_variations().select(  # doctest: +SKIP
            ...     "variation_id", "start_position", "end_position"
            ... ).columns
            ['variation_id', 'start_position', 'end_position']
        """
        return self._extract(
            "protein_sequence_variation",
            "r.variation_id AS variation_id, r.start_position AS start_position, "
            "r.end_position AS end_position, r.note AS note",
            "JOIN protein_sequence_variation r USING (primary_accession)",
            order_by="r.start_position, r.end_position, r.variation_id",
        )

    def extract_isoform_variations(self) -> pl.DataFrame:
        """Return ordered isoform-to-VSP relationships.

        Examples:
            >>> selection.extract_isoform_variations().select(  # doctest: +SKIP
            ...     "isoform_id", "variation_id", "variation_order"
            ... ).columns
            ['isoform_id', 'variation_id', 'variation_order']
        """
        return self._extract(
            "protein_isoform_variation",
            "r.isoform_id AS isoform_id, r.variation_id AS variation_id, "
            "r.variation_order AS variation_order",
            "JOIN protein_isoform_variation r USING (primary_accession)",
            order_by="r.isoform_id, r.variation_order, r.variation_id",
        )

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Return requested identifiers that matched no protein after filtering.

        Examples:
            >>> selection.extract_unmatched_ids().select(  # doctest: +SKIP
            ...     "input_id", "input_namespace", "reason"
            ... ).columns
            ['input_id', 'input_namespace', 'reason']
        """
        frame = self._identifier_matches()
        matched: set[str] = (
            set(frame["input_id"].cast(pl.String).to_list()) if frame.height else set()
        )
        rows = [
            {
                "group_id": group,
                "input_id": input_id,
                "input_namespace": self.namespace,
                "reason": "not_found",
            }
            for group, input_id in self.group_membership
            if input_id not in matched
        ]
        frame = pl.DataFrame(
            rows,
            schema={
                "group_id": pl.String,
                "input_id": pl.String,
                "input_namespace": pl.String,
                "reason": pl.String,
            },
        )
        return frame if self._is_grouped else frame.drop("group_id")

    def _extract(
        self,
        table: str,
        columns: str,
        join: str,
        *,
        condition: str = "",
        parameters: list[object] | None = None,
        order_by: str,
    ) -> pl.DataFrame:
        del table
        matches = self._matches()
        with self.database.connect() as connection:
            connection.execute(
                "CREATE TEMP TABLE _selection("
                "group_id VARCHAR, input_id VARCHAR, input_namespace VARCHAR, "
                "primary_accession VARCHAR)"
            )
            if matches.height:
                connection.executemany(
                    "INSERT INTO _selection VALUES (?, ?, ?, ?)", matches.rows()
                )
            cursor = connection.execute(
                "SELECT DISTINCT s.group_id, s.input_id, s.input_namespace, "
                "s.primary_accession AS primary_accession, "
                f"{columns} FROM _selection s {join}{condition} "
                "ORDER BY s.group_id, s.input_id, s.primary_accession, "
                f"{order_by}",
                parameters or [],
            )
            frame = _cursor_frame(cursor)
        return frame if self._is_grouped else frame.drop("group_id")

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
        cached = self._matches_cache
        if cached is not None:
            return cached
        matches = self._query_identifier_matches()
        object.__setattr__(self, "_matches_cache", matches)
        return matches

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
        if metadata.get("bioextract.metadata_schema_version") != "1":
            raise ValueError("Unsupported UniProt metadata schema version")
        required_v1 = {
            "bioextract.resource_name",
            "bioextract.resource_schema_version",
            "bioextract.source_schema_profile",
            "bioextract.package_version",
            "bioextract.generated_at",
            "bioextract.validation_status",
            "bioextract.validation_issue_count",
            "bioextract.sources",
        }
        missing_v1 = sorted(required_v1 - set(metadata))
        if missing_v1:
            raise ValueError(f"UniProt metadata v1 is missing keys: {missing_v1}")
        validate_duckdb_metadata_v1(connection, metadata)
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


def _polars_type(type_name: str) -> type[pl.DataType]:
    if type_name == "BOOLEAN":
        return pl.Boolean
    if type_name in {"BIGINT", "INTEGER", "UBIGINT"}:
        return pl.Int64
    return pl.String
