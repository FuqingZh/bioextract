from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._publication import (
    validate_duckdb_metadata_v3,
    validate_duckdb_validation_state,
)

from .constant import SCHEMA_VERSION, SOURCE_SCHEMA_PROFILE, RheaNamespace

if TYPE_CHECKING:
    from .rhea import RheaDatabase


_METADATA_TABLES = {
    "metadata",
    "source_file",
    "table_info",
    "column_mapping",
}
_METADATA_SCHEMA_VERSIONS = {"1", "2", "3"}
_EXTERNAL_DATABASE_BY_NAMESPACE: Mapping[RheaNamespace, str] = {
    "ec": "EC",
    "go": "GO",
    "ecocyc": "ECOCYC",
    "kegg_reaction": "KEGG_REACTION",
    "macie": "MACIE",
    "metacyc": "METACYC",
    "reactome": "REACTOME",
}
_NAMESPACE_CAPABILITIES: Mapping[RheaNamespace, frozenset[str]] = {
    "rhea": frozenset({"reaction"}),
    "chebi": frozenset({"reaction", "reaction_participant", "compound"}),
    "uniprot": frozenset({"reaction", "reaction_uniprot"}),
    **{
        namespace: frozenset({"reaction", "reaction_xref"})
        for namespace in _EXTERNAL_DATABASE_BY_NAMESPACE
    },
}
_RHEA_PATTERN = re.compile(r"^(?:RHEA:)?([0-9]+)$", re.IGNORECASE)
_CHEBI_PATTERN = re.compile(r"^(?:CHEBI:)?([0-9]+)$", re.IGNORECASE)
_GO_PATTERN = re.compile(r"^GO:([0-9]{7})$", re.IGNORECASE)

_SCHEMA_MATCH: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
    "InputNamespace": pl.String,
    "RheaId": pl.Int64,
    "MasterId": pl.Int64,
    "Direction": pl.String,
}
_SCHEMA_REACTION: SchemaDict = {
    **_SCHEMA_MATCH,
    "Accession": pl.String,
    "Equation": pl.String,
    "EquationHtml": pl.String,
    "Status": pl.String,
    "IsBalanced": pl.Boolean,
    "IsTransport": pl.Boolean,
    "Comment": pl.String,
    "IsObsolete": pl.Boolean,
    "ReactionSmiles": pl.String,
}
_SCHEMA_PARTICIPANT: SchemaDict = {
    **_SCHEMA_MATCH,
    "ParticipantId": pl.String,
    "CompoundId": pl.String,
    "Side": pl.String,
    "DirectionalRole": pl.String,
    "CoefficientText": pl.String,
    "CoefficientNumeric": pl.Float64,
    "Location": pl.String,
    "RheaCompoundId": pl.Int64,
    "CompoundAccession": pl.String,
    "CompoundType": pl.String,
    "CompoundName": pl.String,
    "Formula": pl.String,
    "ChargeText": pl.String,
    "ChargeNumeric": pl.Int64,
    "ChEBIId": pl.String,
    "UnderlyingChEBIId": pl.String,
    "PolymerizationIndex": pl.String,
}
_SCHEMA_CROSS_REFERENCE: SchemaDict = {
    **_SCHEMA_MATCH,
    "ReferenceDatabase": pl.String,
    "ReferenceId": pl.String,
    "UniProtSection": pl.String,
}
_SCHEMA_PUBLICATION: SchemaDict = {
    **_SCHEMA_MATCH,
    "PubMedId": pl.String,
}
_SCHEMA_RELATIONSHIP: SchemaDict = {
    **_SCHEMA_MATCH,
    "FromReactionId": pl.Int64,
    "ToReactionId": pl.Int64,
    "RelationType": pl.String,
}
_SCHEMA_UNMATCHED: SchemaDict = {
    "GroupId": pl.String,
    "InputId": pl.String,
    "InputNamespace": pl.String,
}


class RheaCapabilityError(RuntimeError):
    """Raised when a partial Rhea publication lacks a required relation."""


@dataclass(frozen=True, slots=True)
class _RheaPublication:
    path: Path
    tables: frozenset[str]
    views: frozenset[str]
    metadata: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class _InputRow:
    group_id: str | None
    input_id: str
    lookup_value: str


@dataclass(frozen=True, slots=True)
class RheaReactionSelection:
    """Deferred Rhea reaction selection with eager extraction terminals.

    The selection only stores normalized input and query policy. DuckDB is
    opened read-only when an ``extract_*`` terminal is called.

    Examples:
        Inspect the declared namespace without executing the query:

        >>> selection.namespace  # doctest: +SKIP
        'chebi'
    """

    database: RheaDatabase = field(repr=False)
    namespace: RheaNamespace
    include_obsolete: bool
    _input_rows: tuple[_InputRow, ...] = field(repr=False)
    _is_grouped: bool = field(repr=False)
    _matches_cache: pl.DataFrame | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )

    def extract_matches(self) -> pl.DataFrame:
        """Resolve input identifiers to exact Rhea reactions.

        Examples:
            Materialize the identifier-to-reaction mapping:

            >>> selection.extract_matches()["RheaId"].head(1).to_list()  # doctest: +SKIP
            [10000]
        """
        frame = self._matches_cache
        if frame is None:
            frame = self._query_matches()
            object.__setattr__(self, "_matches_cache", frame)
        return frame.clone()

    def extract_reactions(self) -> pl.DataFrame:
        """Extract selected reaction facts, including optional reaction SMILES.

        Examples:
            Read exact direction alongside the reaction equation:

            >>> selection.extract_reactions().select(  # doctest: +SKIP
            ...     "RheaId", "Direction", "Equation"
            ... ).columns
            ['RheaId', 'Direction', 'Equation']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(publication, {"reaction"}, operation="extract reactions")
        has_smiles = "reaction_smiles" in publication.tables
        smiles_expression = (
            "smiles.reaction_smiles" if has_smiles else "CAST(NULL AS VARCHAR)"
        )
        smiles_join = (
            "LEFT JOIN reaction_smiles AS smiles USING (rhea_id)" if has_smiles else ""
        )
        query = f"""
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                reaction.rhea_id AS "RheaId",
                reaction.master_id AS "MasterId",
                reaction.direction AS "Direction",
                reaction.accession AS "Accession",
                reaction.equation AS "Equation",
                reaction.equation_html AS "EquationHtml",
                reaction.status AS "Status",
                reaction.is_balanced AS "IsBalanced",
                reaction.is_transport AS "IsTransport",
                reaction.comment AS "Comment",
                reaction.is_obsolete AS "IsObsolete",
                {smiles_expression} AS "ReactionSmiles"
            FROM _selected_reaction AS selected
            JOIN reaction USING (rhea_id)
            {smiles_join}
            ORDER BY "GroupId", "InputId", "RheaId"
        """
        return self._query_selected(query, schema=_SCHEMA_REACTION)

    def extract_participants(self) -> pl.DataFrame:
        """Extract participants with compound facts and direction-aware roles.

        Examples:
            Read retained sides and nullable direction-specific roles:

            >>> selection.extract_participants().select(  # doctest: +SKIP
            ...     "Side", "DirectionalRole", "ChEBIId"
            ... ).columns
            ['Side', 'DirectionalRole', 'ChEBIId']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_participant", "compound"},
            operation="extract reaction participants",
        )
        if "reaction_participant_direction" not in publication.views:
            raise RheaCapabilityError(
                "Rhea publication lacks view required to extract reaction "
                "participants: reaction_participant_direction"
            )
        query = """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                participant.rhea_id AS "RheaId",
                participant.master_id AS "MasterId",
                participant.direction AS "Direction",
                participant.participant_id AS "ParticipantId",
                participant.compound_id AS "CompoundId",
                participant.side AS "Side",
                participant.directional_role AS "DirectionalRole",
                participant.coefficient_text AS "CoefficientText",
                participant.coefficient_numeric AS "CoefficientNumeric",
                participant.location AS "Location",
                compound.rhea_compound_id AS "RheaCompoundId",
                compound.public_accession AS "CompoundAccession",
                compound.compound_type AS "CompoundType",
                compound.name AS "CompoundName",
                compound.formula AS "Formula",
                compound.charge_text AS "ChargeText",
                compound.charge_numeric AS "ChargeNumeric",
                compound.chebi_id AS "ChEBIId",
                compound.underlying_chebi_id AS "UnderlyingChEBIId",
                compound.polymerization_index AS "PolymerizationIndex"
            FROM _selected_reaction AS selected
            JOIN reaction_participant_direction AS participant USING (rhea_id)
            LEFT JOIN compound USING (compound_id)
            ORDER BY
                "GroupId", "InputId", "RheaId", "Side", "ParticipantId"
        """
        return self._query_selected(query, schema=_SCHEMA_PARTICIPANT)

    def extract_cross_references(self) -> pl.DataFrame:
        """Extract Rhea-owned external and UniProt reaction mappings.

        Examples:
            Materialize the normalized reference relation:

            >>> selection.extract_cross_references().select(  # doctest: +SKIP
            ...     "ReferenceDatabase", "ReferenceId"
            ... ).columns
            ['ReferenceDatabase', 'ReferenceId']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        has_xref = "reaction_xref" in publication.tables
        has_uniprot = "reaction_uniprot" in publication.tables
        if not has_xref and not has_uniprot:
            raise RheaCapabilityError(
                "Rhea publication lacks relations required to extract "
                "cross-references: reaction_xref or reaction_uniprot"
            )
        branches: list[str] = []
        if has_xref:
            branches.append(
                """
                SELECT
                    selected.group_id,
                    selected.input_id,
                    selected.input_namespace,
                    selected.rhea_id,
                    selected.master_id,
                    selected.direction,
                    xref.database_name AS reference_database,
                    xref.external_id AS reference_id,
                    CAST(NULL AS VARCHAR) AS uniprot_section
                FROM _selected_reaction AS selected
                JOIN reaction_xref AS xref USING (rhea_id)
                """
            )
        if has_uniprot:
            branches.append(
                """
                SELECT
                    selected.group_id,
                    selected.input_id,
                    selected.input_namespace,
                    selected.rhea_id,
                    selected.master_id,
                    selected.direction,
                    'UniProt' AS reference_database,
                    mapping.uniprot_id AS reference_id,
                    mapping.uniprot_section
                FROM _selected_reaction AS selected
                JOIN reaction_uniprot AS mapping USING (rhea_id)
                """
            )
        query = f"""
            SELECT
                group_id AS "GroupId",
                input_id AS "InputId",
                input_namespace AS "InputNamespace",
                rhea_id AS "RheaId",
                master_id AS "MasterId",
                direction AS "Direction",
                reference_database AS "ReferenceDatabase",
                reference_id AS "ReferenceId",
                uniprot_section AS "UniProtSection"
            FROM ({" UNION ALL ".join(branches)})
            ORDER BY
                "GroupId", "InputId", "RheaId",
                "ReferenceDatabase", "ReferenceId"
        """
        return self._query_selected(query, schema=_SCHEMA_CROSS_REFERENCE)

    def extract_publications(self) -> pl.DataFrame:
        """Extract PubMed references owned by selected Rhea reactions.

        Examples:
            Read PubMed IDs while retaining reaction lineage:

            >>> selection.extract_publications().select(  # doctest: +SKIP
            ...     "RheaId", "PubMedId"
            ... ).columns
            ['RheaId', 'PubMedId']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_publication"},
            operation="extract reaction publications",
        )
        query = """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.rhea_id AS "RheaId",
                selected.master_id AS "MasterId",
                selected.direction AS "Direction",
                publication.pubmed_id AS "PubMedId"
            FROM _selected_reaction AS selected
            JOIN reaction_publication AS publication USING (rhea_id)
            ORDER BY "GroupId", "InputId", "RheaId", "PubMedId"
        """
        return self._query_selected(query, schema=_SCHEMA_PUBLICATION)

    def extract_relationships(self) -> pl.DataFrame:
        """Extract hierarchy edges touching selected master reactions.

        Examples:
            Preserve each directed Rhea hierarchy edge:

            >>> selection.extract_relationships().select(  # doctest: +SKIP
            ...     "FromReactionId", "ToReactionId", "RelationType"
            ... ).columns
            ['FromReactionId', 'ToReactionId', 'RelationType']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_relationship"},
            operation="extract reaction relationships",
        )
        query = """
            SELECT DISTINCT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.rhea_id AS "RheaId",
                selected.master_id AS "MasterId",
                selected.direction AS "Direction",
                relation.from_reaction_id AS "FromReactionId",
                relation.to_reaction_id AS "ToReactionId",
                relation.relation_type AS "RelationType"
            FROM _selected_reaction AS selected
            JOIN reaction_relationship AS relation
              ON relation.from_reaction_id = selected.master_id
              OR relation.to_reaction_id = selected.master_id
            ORDER BY
                "GroupId", "InputId", "RheaId",
                "FromReactionId", "ToReactionId", "RelationType"
        """
        return self._query_selected(query, schema=_SCHEMA_RELATIONSHIP)

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Return normalized input identifiers with no matching Rhea reaction.

        Examples:
            Audit identifiers that did not resolve:

            >>> selection.extract_unmatched_ids().columns  # doctest: +SKIP
            ['InputId', 'InputNamespace']
        """
        matched = self.extract_matches()
        matched_keys = {
            (row.get("GroupId"), row["InputId"]) for row in matched.to_dicts()
        }
        rows = [
            (row.group_id, row.input_id, self.namespace)
            for row in self._input_rows
            if (row.group_id, row.input_id) not in matched_keys
        ]
        return self._finalize(
            pl.DataFrame(rows, schema=_SCHEMA_UNMATCHED, orient="row"),
            schema=_SCHEMA_UNMATCHED,
        )

    def _query_matches(self) -> pl.DataFrame:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            set(_NAMESPACE_CAPABILITIES[self.namespace]),
            operation=f"select reactions by {self.namespace}",
        )
        obsolete_filter = (
            "" if self.include_obsolete else "AND NOT reaction.is_obsolete"
        )
        if self.namespace == "rhea":
            join = (
                "JOIN reaction "
                "ON CAST(reaction.rhea_id AS VARCHAR) = input.lookup_value"
            )
        elif self.namespace == "chebi":
            join = """
                JOIN compound
                  ON compound.chebi_id = input.lookup_value
                JOIN reaction_participant AS participant
                  ON participant.compound_id = compound.compound_id
                JOIN reaction
                  ON reaction.master_id = participant.master_id
            """
        elif self.namespace == "uniprot":
            join = """
                JOIN reaction_uniprot AS mapping
                  ON mapping.uniprot_id = input.lookup_value
                JOIN reaction USING (rhea_id)
            """
        else:
            database_name = _EXTERNAL_DATABASE_BY_NAMESPACE[self.namespace]
            join = f"""
                JOIN reaction_xref AS mapping
                  ON mapping.external_id = input.lookup_value
                 AND mapping.database_name = '{database_name}'
                JOIN reaction USING (rhea_id)
            """
        query = f"""
            SELECT DISTINCT
                input.group_id AS "GroupId",
                input.input_id AS "InputId",
                '{self.namespace}' AS "InputNamespace",
                reaction.rhea_id AS "RheaId",
                reaction.master_id AS "MasterId",
                reaction.direction AS "Direction"
            FROM _input_id AS input
            {join}
            WHERE true {obsolete_filter}
            ORDER BY "GroupId", "InputId", "RheaId"
        """
        with _connect(publication) as connection:
            _create_input_table(connection, self._input_rows)
            frame = _fetch_frame(connection, query, schema=_SCHEMA_MATCH)
        return self._finalize(frame, schema=_SCHEMA_MATCH)

    def _query_selected(
        self,
        query: str,
        *,
        schema: SchemaDict,
    ) -> pl.DataFrame:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        matches = self.extract_matches()
        with _connect(publication) as connection:
            _create_selected_table(connection, matches, grouped=self._is_grouped)
            frame = _fetch_frame(connection, query, schema=schema)
        return self._finalize(frame, schema=schema)

    def _finalize(
        self,
        frame: pl.DataFrame,
        *,
        schema: SchemaDict,
    ) -> pl.DataFrame:
        expected = list(schema)
        if not self._is_grouped:
            expected.remove("GroupId")
            frame = frame.drop("GroupId")
        return frame.select(expected)


def open_rhea_publication(path: str | Path) -> _RheaPublication:
    """Validate and describe one bioextract Rhea DuckDB publication."""
    publication_path = Path(path)
    if not publication_path.exists():
        raise FileNotFoundError(publication_path)
    if not publication_path.is_file():
        raise ValueError(f"Rhea DuckDB path is not a file: {publication_path}")

    try:
        connection = duckdb.connect(str(publication_path), read_only=True)
    except duckdb.Error as error:
        raise ValueError(
            f"Cannot open Rhea DuckDB publication: {publication_path}"
        ) from error
    try:
        metadata_tables = {
            str(row[0])
            for row in connection.execute(
                """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = '_bioextract'
                  AND table_type = 'BASE TABLE'
                """
            ).fetchall()
        }
        missing_metadata_tables = sorted(_METADATA_TABLES - metadata_tables)
        if missing_metadata_tables:
            raise ValueError(
                "DuckDB file is missing bioextract metadata tables: "
                f"{missing_metadata_tables}"
            )

        metadata = {
            str(key): str(value)
            for key, value in connection.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
        }
        metadata_schema_version = metadata.get("bioextract.metadata_schema_version")
        if metadata_schema_version not in _METADATA_SCHEMA_VERSIONS:
            raise ValueError(
                "Unsupported Rhea publication metadata: "
                f"bioextract.metadata_schema_version={metadata_schema_version!r}"
            )
        if metadata_schema_version == "3":
            validate_duckdb_metadata_v3(connection, metadata)
            required_v3 = {
                "bioextract.resource_name",
                "bioextract.resource_schema_version",
                "bioextract.source_schema_profile",
                "bioextract.package_version",
                "bioextract.generated_at",
                "bioextract.validation_status",
                "bioextract.validation_issue_count",
                "bioextract.sources",
            }
            missing_v3 = sorted(required_v3 - set(metadata))
            if missing_v3:
                raise ValueError(f"Rhea metadata v3 is missing keys: {missing_v3}")
            if (
                metadata.get("bioextract.source_schema_profile")
                != SOURCE_SCHEMA_PROFILE
            ):
                raise ValueError("Unsupported Rhea source schema profile")
        if metadata_schema_version == "2":
            validate_duckdb_validation_state(connection, metadata)
        resource_schema_key = (
            "bioextract.resource_schema_version"
            if metadata_schema_version == "3"
            else "bioextract.schema_version"
        )
        expected_metadata = {
            "bioextract.resource_name": "rhea",
            resource_schema_key: SCHEMA_VERSION,
        }
        for key, expected in expected_metadata.items():
            observed = metadata.get(key)
            if observed != expected:
                raise ValueError(
                    "Unsupported Rhea publication metadata: "
                    f"{key}={observed!r}, expected {expected!r}"
                )

        tables = frozenset(
            str(row[0])
            for row in connection.execute(
                """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = 'main'
                  AND table_type = 'BASE TABLE'
                """
            ).fetchall()
        )
        recorded_counts = {
            str(name): int(count)
            for name, count in connection.execute(
                "SELECT table_name, row_count FROM _bioextract.table_info"
            ).fetchall()
        }
        if set(recorded_counts) != set(tables):
            raise ValueError(
                "Rhea publication table inventory does not match main tables"
            )
        for table_name, expected_count in recorded_counts.items():
            observed_row = connection.execute(
                f'SELECT count(*) FROM "{table_name}"'
            ).fetchone()
            observed_count = 0 if observed_row is None else int(observed_row[0])
            if observed_count != expected_count:
                raise ValueError(
                    "Rhea publication row count does not match metadata: "
                    f"{table_name}={observed_count}, expected {expected_count}"
                )
        if "compound" in tables:
            chebi_type_row = connection.execute(
                """
                SELECT data_type
                FROM information_schema.columns
                WHERE table_schema = 'main'
                  AND table_name = 'compound'
                  AND column_name = 'chebi_id'
                """
            ).fetchone()
            if chebi_type_row is None or str(chebi_type_row[0]) != "VARCHAR":
                raise ValueError(
                    "Unsupported Rhea v1 physical schema: "
                    "compound.chebi_id must be VARCHAR CURIE"
                )
            malformed_row = connection.execute(
                """
                SELECT count(*)
                FROM compound
                WHERE chebi_id IS NOT NULL
                  AND NOT regexp_full_match(chebi_id, 'CHEBI:[0-9]+')
                """
            ).fetchone()
            if malformed_row is None:
                raise ValueError("Cannot validate Rhea ChEBI identifiers")
            malformed_count = int(malformed_row[0])
            if malformed_count:
                raise ValueError(
                    "Unsupported Rhea v1 physical schema: "
                    "compound.chebi_id contains non-CURIE values"
                )
        views = frozenset(
            str(row[0])
            for row in connection.execute(
                """
                SELECT table_name
                FROM information_schema.views
                WHERE table_schema = 'main'
                """
            ).fetchall()
        )
    except duckdb.Error as error:
        raise ValueError(
            f"Invalid bioextract Rhea publication: {publication_path}"
        ) from error
    finally:
        connection.close()
    return _RheaPublication(
        path=publication_path,
        tables=tables,
        views=views,
        metadata=metadata,
    )


def create_selection(
    database: RheaDatabase,
    ids: Iterable[str],
    *,
    namespace: RheaNamespace,
    include_obsolete: bool,
) -> RheaReactionSelection:
    publication = database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
    normalized_namespace = _validate_namespace(namespace)
    _require_tables(
        publication,
        set(_NAMESPACE_CAPABILITIES[normalized_namespace]),
        operation=f"select reactions by {normalized_namespace}",
    )
    rows = tuple(
        _InputRow(group_id=None, input_id=input_id, lookup_value=lookup_value)
        for input_id, lookup_value in _normalize_ids(ids, normalized_namespace)
    )
    return RheaReactionSelection(
        database=database,
        namespace=normalized_namespace,
        include_obsolete=bool(include_obsolete),
        _input_rows=rows,
        _is_grouped=False,
    )


def create_group_selection(
    database: RheaDatabase,
    ids_by_group: Mapping[str, Iterable[str]],
    *,
    namespace: RheaNamespace,
    include_obsolete: bool,
) -> RheaReactionSelection:
    publication = database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
    normalized_namespace = _validate_namespace(namespace)
    _require_tables(
        publication,
        set(_NAMESPACE_CAPABILITIES[normalized_namespace]),
        operation=f"select reactions by {normalized_namespace}",
    )
    rows: list[_InputRow] = []
    for raw_group_id, ids in ids_by_group.items():
        group_id = str(raw_group_id).strip()
        if not group_id:
            raise ValueError("Rhea group IDs must be non-empty after normalization")
        rows.extend(
            _InputRow(
                group_id=group_id,
                input_id=input_id,
                lookup_value=lookup_value,
            )
            for input_id, lookup_value in _normalize_ids(ids, normalized_namespace)
        )
    rows = sorted(
        set(rows),
        key=lambda row: (row.group_id or "", row.input_id, row.lookup_value),
    )
    return RheaReactionSelection(
        database=database,
        namespace=normalized_namespace,
        include_obsolete=bool(include_obsolete),
        _input_rows=tuple(rows),
        _is_grouped=True,
    )


def _validate_namespace(namespace: str) -> RheaNamespace:
    if namespace not in _NAMESPACE_CAPABILITIES:
        supported = ", ".join(_NAMESPACE_CAPABILITIES)
        raise ValueError(
            f"Unsupported Rhea namespace: {namespace!r}; expected one of {supported}"
        )
    return namespace  # type: ignore[return-value]


def _normalize_ids(
    ids: Iterable[str],
    namespace: RheaNamespace,
) -> list[tuple[str, str]]:
    normalized: set[tuple[str, str]] = set()
    for value in ids:
        text = str(value).strip()
        if not text:
            continue
        if namespace == "rhea":
            match = _RHEA_PATTERN.fullmatch(text)
            if match is None:
                raise ValueError(f"Invalid Rhea identifier: {value!r}")
            lookup = str(int(match.group(1)))
            normalized.add((f"RHEA:{lookup}", lookup))
        elif namespace == "chebi":
            match = _CHEBI_PATTERN.fullmatch(text)
            if match is None:
                raise ValueError(f"Invalid ChEBI identifier: {value!r}")
            lookup = f"CHEBI:{int(match.group(1))}"
            normalized.add((lookup, lookup))
        elif namespace == "go":
            match = _GO_PATTERN.fullmatch(text)
            if match is None:
                raise ValueError(f"Invalid GO identifier: {value!r}")
            input_id = f"GO:{match.group(1)}"
            normalized.add((input_id, input_id))
        elif namespace == "uniprot" and "|" in text:
            parts = text.split("|")
            if len(parts) >= 2 and parts[1].strip():
                accession = parts[1].strip()
                normalized.add((accession, accession))
            else:
                raise ValueError(f"Invalid UniProt identifier: {value!r}")
        else:
            normalized.add((text, text))
    return sorted(normalized)


def _require_tables(
    publication: _RheaPublication,
    required: set[str],
    *,
    operation: str,
) -> None:
    missing = sorted(required - set(publication.tables))
    if missing:
        raise RheaCapabilityError(
            f"Rhea publication cannot {operation}; missing relations: {missing}"
        )


def _connect(
    publication: _RheaPublication,
) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(str(publication.path), read_only=True)


def _create_input_table(
    connection: duckdb.DuckDBPyConnection,
    rows: tuple[_InputRow, ...],
) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _input_id (
            group_id VARCHAR,
            input_id VARCHAR NOT NULL,
            lookup_value VARCHAR NOT NULL
        )
        """
    )
    if rows:
        connection.executemany(
            "INSERT INTO _input_id VALUES (?, ?, ?)",
            [(row.group_id, row.input_id, row.lookup_value) for row in rows],
        )


def _create_selected_table(
    connection: duckdb.DuckDBPyConnection,
    matches: pl.DataFrame,
    *,
    grouped: bool,
) -> None:
    connection.execute(
        """
        CREATE TEMP TABLE _selected_reaction (
            group_id VARCHAR,
            input_id VARCHAR NOT NULL,
            input_namespace VARCHAR NOT NULL,
            rhea_id BIGINT NOT NULL,
            master_id BIGINT,
            direction VARCHAR
        )
        """
    )
    rows = [
        (
            row["GroupId"] if grouped else None,
            row["InputId"],
            row["InputNamespace"],
            row["RheaId"],
            row["MasterId"],
            row["Direction"],
        )
        for row in matches.to_dicts()
    ]
    if rows:
        connection.executemany(
            "INSERT INTO _selected_reaction VALUES (?, ?, ?, ?, ?, ?)",
            rows,
        )


def _fetch_frame(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    *,
    schema: SchemaDict,
) -> pl.DataFrame:
    rows = connection.execute(query).fetchall()
    return pl.DataFrame(rows, schema=schema, orient="row", strict=False)
