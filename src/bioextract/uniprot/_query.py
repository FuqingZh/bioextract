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
        ...     "UniProtId", "Accession", "IsPrimaryAccession"
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
            ...     "UniProtId", "EntryName", "TaxonId"
            ... ).columns
            ['UniProtId', 'EntryName', 'TaxonId']
        """
        return self._extract(
            "protein",
            "p.entry_name AS EntryName, p.is_reviewed AS IsReviewed, "
            "p.taxon_id AS TaxonId, p.protein_existence AS ProteinExistence, "
            "p.sequence_length AS SequenceLength, "
            "p.molecular_weight AS MolecularWeight, "
            "p.sequence_version AS SequenceVersion, "
            "p.entry_version AS EntryVersion",
            "JOIN protein p USING (primary_accession)",
            order_by="p.primary_accession",
        )

    def extract_accessions(self) -> pl.DataFrame:
        """Return ordered primary and secondary accessions for every match.

        Examples:
            >>> selection.extract_accessions().select(  # doctest: +SKIP
            ...     "Accession", "IsPrimaryAccession"
            ... ).columns
            ['Accession', 'IsPrimaryAccession']
        """
        return self._extract(
            "protein_accession",
            "r.accession AS Accession, r.accession_order AS AccessionOrder, "
            "r.is_primary AS IsPrimaryAccession",
            "JOIN protein_accession r USING (primary_accession)",
            order_by="r.accession_order, r.accession",
        )

    def extract_protein_names(self) -> pl.DataFrame:
        """Return official protein names in source-defined name order.

        Examples:
            >>> selection.extract_protein_names().select(  # doctest: +SKIP
            ...     "NameType", "ProteinName"
            ... ).columns
            ['NameType', 'ProteinName']
        """
        return self._extract(
            "protein_name",
            "r.name_type AS NameType, r.name AS ProteinName, r.name_order AS NameOrder",
            "JOIN protein_name r USING (primary_accession)",
            order_by="r.name_order, r.name",
        )

    def extract_gene_names(self) -> pl.DataFrame:
        """Return parsed gene names in source-defined name order.

        Examples:
            >>> selection.extract_gene_names().select(  # doctest: +SKIP
            ...     "NameType", "GeneName"
            ... ).columns
            ['NameType', 'GeneName']
        """
        return self._extract(
            "gene_name",
            "r.name_type AS NameType, r.name AS GeneName, r.name_order AS NameOrder",
            "JOIN gene_name r USING (primary_accession)",
            order_by="r.name_order, r.name",
        )

    def extract_ec_numbers(self) -> pl.DataFrame:
        """Return distinct EC annotations ordered by EC number.

        Examples:
            >>> selection.extract_ec_numbers().select("EcNumber").columns  # doctest: +SKIP
            ['EcNumber']
        """
        return self._extract(
            "protein_ec_number",
            "r.ec_number AS ECNumber",
            "JOIN protein_ec_number r USING (primary_accession)",
            order_by="r.ec_number",
        )

    def extract_go_annotations(self) -> pl.DataFrame:
        """Return GO annotations ordered by GO identifier and aspect.

        Examples:
            >>> selection.extract_go_annotations().select(  # doctest: +SKIP
            ...     "GOId", "Aspect", "EvidenceCode"
            ... ).columns
            ['GOId', 'Aspect', 'EvidenceCode']
        """
        return self._extract(
            "protein_go_annotation",
            "r.go_id AS GOId, r.aspect AS Aspect, r.term_name AS TermName, "
            "r.evidence_code AS EvidenceCode, r.evidence_source AS EvidenceSource",
            "JOIN protein_go_annotation r USING (primary_accession)",
            order_by="r.go_id, r.aspect",
        )

    def extract_cross_references(
        self, databases: Iterable[str] | None = None
    ) -> pl.DataFrame:
        """Return cross-references, optionally limited to database names.

        Examples:
            >>> selection.extract_cross_references(["PDB"]).select(  # doctest: +SKIP
            ...     "Database", "ExternalId"
            ... ).columns
            ['Database', 'ExternalId']
        """
        values = tuple(dict.fromkeys(databases or ()))
        condition = ""
        parameters: list[object] = []
        if values:
            condition = f" WHERE r.database IN ({','.join('?' for _ in values)})"
            parameters.extend(values)
        return self._extract(
            "protein_cross_reference",
            "r.database AS Database, r.external_id AS ExternalId, "
            "r.properties AS Properties, r.isoform_id AS IsoformId",
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
            ...     "CommentType", "CommentText"
            ... ).columns
            ['CommentType', 'CommentText']
        """
        values = tuple(dict.fromkeys(comment_types or ()))
        condition = ""
        parameters: list[object] = []
        if values:
            condition = f" WHERE r.comment_type IN ({','.join('?' for _ in values)})"
            parameters.extend(values)
        return self._extract(
            "protein_comment",
            "r.comment_id AS CommentId, r.comment_type AS CommentType, "
            "r.comment_text AS CommentText",
            "JOIN protein_comment r USING (primary_accession)",
            condition=condition,
            parameters=parameters,
            order_by="r.comment_id",
        )

    def extract_subcellular_locations(self) -> pl.DataFrame:
        """Return individually parsed locations and their optional notes.

        Examples:
            >>> selection.extract_subcellular_locations().select(  # doctest: +SKIP
            ...     "SubcellularLocation", "SubcellularLocationNote"
            ... ).columns
            ['SubcellularLocation', 'SubcellularLocationNote']
        """
        return self._extract(
            "protein_subcellular_location",
            "r.location AS SubcellularLocation, r.note AS SubcellularLocationNote",
            "JOIN protein_subcellular_location r USING (primary_accession)",
            order_by="r.comment_id, r.location",
        )

    def extract_keywords(self) -> pl.DataFrame:
        """Return keywords in source-defined keyword order.

        Examples:
            >>> selection.extract_keywords().select(  # doctest: +SKIP
            ...     "Keyword", "KeywordOrder"
            ... ).columns
            ['Keyword', 'KeywordOrder']
        """
        return self._extract(
            "protein_keyword",
            "r.keyword AS Keyword, r.keyword_order AS KeywordOrder",
            "JOIN protein_keyword r USING (primary_accession)",
            order_by="r.keyword_order, r.keyword",
        )

    def extract_sequences(self, sequence_type: str = "canonical") -> pl.DataFrame:
        """Return canonical, isoform, or all sequences in stable type order.

        Examples:
            >>> selection.extract_sequences("all").select(  # doctest: +SKIP
            ...     "SequenceId", "SequenceType", "SHA256"
            ... ).columns
            ['SequenceId', 'SequenceType', 'SHA256']
        """
        if sequence_type not in {"canonical", "isoform", "all"}:
            raise ValueError("sequence_type must be canonical, isoform, or all")
        condition = "" if sequence_type == "all" else " WHERE r.sequence_type = ?"
        parameters: list[object] = [] if sequence_type == "all" else [sequence_type]
        return self._extract(
            "protein_sequence",
            "r.sequence_id AS SequenceId, r.sequence_type AS SequenceType, "
            "r.sequence AS Sequence, r.length AS Length, r.crc64 AS CRC64, "
            "r.sha256 AS SHA256",
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
            ...     "IsoformId", "SequenceStatus", "SequenceId"
            ... ).columns
            ['IsoformId', 'SequenceStatus', 'SequenceId']
        """
        return self._extract(
            "protein_isoform",
            "r.isoform_id AS IsoformId, r.name AS IsoformName, "
            "r.isoform_order AS IsoformOrder, "
            "r.sequence_status AS SequenceStatus, r.sequence_id AS SequenceId",
            "JOIN protein_isoform r USING (primary_accession)",
            order_by="r.isoform_order, r.isoform_id",
        )

    def extract_isoform_identifiers(self) -> pl.DataFrame:
        """Return every current or old official IsoId for each product.

        Examples:
            >>> selection.extract_isoform_identifiers().select(  # doctest: +SKIP
            ...     "IsoformId", "Identifier", "IsMain"
            ... ).columns
            ['IsoformId', 'Identifier', 'IsMain']
        """
        return self._extract(
            "protein_isoform_identifier",
            "r.isoform_id AS IsoformId, r.identifier AS Identifier, "
            "r.identifier_order AS IdentifierOrder, r.is_main AS IsMain",
            "JOIN protein_isoform_identifier r USING (primary_accession)",
            order_by="r.isoform_id, r.identifier_order",
        )

    def extract_sequence_variations(self) -> pl.DataFrame:
        """Return DAT VAR_SEQ features ordered by coordinates and VSP ID.

        Examples:
            >>> selection.extract_sequence_variations().select(  # doctest: +SKIP
            ...     "VariationId", "StartPosition", "EndPosition"
            ... ).columns
            ['VariationId', 'StartPosition', 'EndPosition']
        """
        return self._extract(
            "protein_sequence_variation",
            "r.variation_id AS VariationId, r.start_position AS StartPosition, "
            "r.end_position AS EndPosition, r.note AS Note",
            "JOIN protein_sequence_variation r USING (primary_accession)",
            order_by="r.start_position, r.end_position, r.variation_id",
        )

    def extract_isoform_variations(self) -> pl.DataFrame:
        """Return ordered isoform-to-VSP relationships.

        Examples:
            >>> selection.extract_isoform_variations().select(  # doctest: +SKIP
            ...     "IsoformId", "VariationId", "VariationOrder"
            ... ).columns
            ['IsoformId', 'VariationId', 'VariationOrder']
        """
        return self._extract(
            "protein_isoform_variation",
            "r.isoform_id AS IsoformId, r.variation_id AS VariationId, "
            "r.variation_order AS VariationOrder",
            "JOIN protein_isoform_variation r USING (primary_accession)",
            order_by="r.isoform_id, r.variation_order, r.variation_id",
        )

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Return requested identifiers that matched no protein after filtering.

        Examples:
            >>> selection.extract_unmatched_ids().select(  # doctest: +SKIP
            ...     "InputId", "InputNamespace", "Reason"
            ... ).columns
            ['InputId', 'InputNamespace', 'Reason']
        """
        frame = self._identifier_matches()
        matched: set[str] = (
            set(frame["InputId"].cast(pl.String).to_list()) if frame.height else set()
        )
        rows = [
            {
                "GroupId": group,
                "InputId": input_id,
                "InputNamespace": self.namespace,
                "Reason": "not_found",
            }
            for group, input_id in self.group_membership
            if input_id not in matched
        ]
        frame = pl.DataFrame(
            rows,
            schema={
                "GroupId": pl.String,
                "InputId": pl.String,
                "InputNamespace": pl.String,
                "Reason": pl.String,
            },
        )
        return frame if self._is_grouped else frame.drop("GroupId")

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
                "GroupId VARCHAR, InputId VARCHAR, InputNamespace VARCHAR, "
                "primary_accession VARCHAR)"
            )
            if matches.height:
                connection.executemany(
                    "INSERT INTO _selection VALUES (?, ?, ?, ?)", matches.rows()
                )
            cursor = connection.execute(
                "SELECT DISTINCT s.GroupId, s.InputId, s.InputNamespace, "
                "s.primary_accession AS UniProtId, "
                f"{columns} FROM _selection s {join}{condition} "
                "ORDER BY s.GroupId, s.InputId, s.primary_accession, "
                f"{order_by}",
                parameters or [],
            )
            frame = _cursor_frame(cursor)
        return frame if self._is_grouped else frame.drop("GroupId")

    def _matches(self) -> pl.DataFrame:
        matches = self._identifier_matches()
        membership = pl.DataFrame(
            self.group_membership,
            schema={"GroupId": pl.String, "InputId": pl.String},
            orient="row",
        )
        if membership.is_empty() or matches.is_empty():
            return pl.DataFrame(
                schema={
                    "GroupId": pl.String,
                    "InputId": pl.String,
                    "InputNamespace": pl.String,
                    "primary_accession": pl.String,
                }
            )
        return (
            membership.join(matches, on="InputId", how="inner")
            .select("GroupId", "InputId", "InputNamespace", "primary_accession")
            .sort("GroupId", "InputId", "primary_accession")
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
            {"InputId": self.input_ids},
            schema={"InputId": pl.String},
        )
        with self.database.connect() as connection:
            connection.execute("CREATE TEMP TABLE _inputs(InputId VARCHAR)")
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
                SELECT DISTINCT i.InputId, ? AS InputNamespace,
                       p.primary_accession
                FROM _inputs i
                JOIN protein_identifier p
                  ON p.identifier = i.InputId AND p.namespace = ?
                JOIN protein pr USING (primary_accession)
                {taxon_filter}
                ORDER BY i.InputId, p.primary_accession
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
            valid_all_taxa = scope == {"all_taxa": True}
            scope_mapping = (
                cast(dict[str, object], scope) if isinstance(scope, dict) else {}
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
                    isinstance(taxon_id, str) and taxon_id for taxon_id in taxon_ids
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
