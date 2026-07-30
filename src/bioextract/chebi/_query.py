from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import duckdb
import polars as pl

from bioextract._publication import (
    validate_duckdb_metadata_v3,
    validate_duckdb_validation_state,
)

if TYPE_CHECKING:
    from .chebi import ChEBIDatabase

_METADATA_VERSIONS = {"1", "2", "3"}
SOURCE_SCHEMA_PROFILE = "chebi-release-bundle-v1"
_BASE_METADATA_TABLES = {
    "metadata",
    "source_file",
    "table_info",
    "column_mapping",
}
_CHEBI_ID = re.compile(r"^(?:CHEBI:)?([0-9]+)$", re.IGNORECASE)
_BUILTIN_NAMESPACES = {"chebi", "inchi", "inchi_key"}
_REASONS = {
    "not_found",
    "below_min_star_rating",
    "obsolete_excluded",
    "invalid_canonical_target",
}


class ChEBICapabilityError(RuntimeError):
    """Raised when a ChEBI publication lacks data required by an operation."""


@dataclass(frozen=True, slots=True)
class _ChEBIPublication:
    path: Path
    tables: frozenset[str]
    metadata: Mapping[str, str]
    namespaces: frozenset[str]


@dataclass(frozen=True, slots=True)
class _InputRow:
    group_id: str | None
    input_id: str
    lookup_value: str


@dataclass(frozen=True, slots=True)
class ChEBICompoundSelection:
    """Deferred identifier selection over a ChEBI publication.

    Examples:
        Inspect policy without executing the publication query:

        >>> selection.namespace  # doctest: +SKIP
        'chebi'
    """

    database: ChEBIDatabase = field(repr=False)
    namespace: str
    min_star_rating: int
    include_obsolete: bool
    _input_rows: tuple[_InputRow, ...] = field(repr=False)
    _is_grouped: bool = field(repr=False)
    _candidate_cache: pl.DataFrame | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def extract_matches(self) -> pl.DataFrame:
        """Return identifier-to-canonical-compound matches.

        Examples:
            Observe the canonical shared identifier:

            >>> selection.extract_matches()["ChEBIId"].head(1).to_list()  # doctest: +SKIP
            ['CHEBI:15377']
        """
        candidates = self._candidates()
        if candidates.is_empty():
            frame = pl.DataFrame(
                schema={
                    "GroupId": pl.String,
                    "InputId": pl.String,
                    "InputNamespace": pl.String,
                    "ChEBIId": pl.String,
                    "MatchType": pl.String,
                    "MatchedValue": pl.String,
                }
            )
        else:
            frame = (
                candidates.filter(
                    pl.col("_target_exists")
                    & (
                        pl.col("_is_obsolete").fill_null(False).not_()
                        | pl.lit(self.include_obsolete)
                    )
                    & (pl.col("_star_rating").fill_null(0) >= self.min_star_rating)
                )
                .select(
                    "GroupId",
                    "InputId",
                    pl.lit(self.namespace).alias("InputNamespace"),
                    pl.col("_chebi_id").alias("ChEBIId"),
                    pl.col("_match_type").alias("MatchType"),
                    pl.col("_matched_value").alias("MatchedValue"),
                )
                .unique()
                .sort(["GroupId", "InputId", "ChEBIId"])
            )
        return self._finalize(frame)

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Return inputs without a policy-accepted match and a stable reason.

        Examples:
            Inspect stable reason codes separately from successful matches:

            >>> unmatched = selection.extract_unmatched_ids()  # doctest: +SKIP
            >>> "Reason" in unmatched.columns  # doctest: +SKIP
            True
        """
        matches = (
            {
                (row.get("GroupId"), row["InputId"])
                for row in self.extract_matches()
                .with_columns(pl.lit(None, dtype=pl.String).alias("GroupId"))
                .to_dicts()
            }
            if not self._is_grouped
            else {
                (row["GroupId"], row["InputId"])
                for row in self.extract_matches().to_dicts()
            }
        )
        candidate_rows = self._candidates().to_dicts()
        by_input: dict[tuple[str | None, str], list[dict[str, object]]] = {}
        for row in candidate_rows:
            by_input.setdefault((row["GroupId"], str(row["InputId"])), []).append(row)
        issue_inputs = self._invalid_target_inputs()
        output: list[dict[str, str | None]] = []
        for row in self._input_rows:
            key = (row.group_id, row.input_id)
            if key in matches:
                continue
            candidates = by_input.get(key, [])
            if key in issue_inputs or any(
                not bool(candidate["_target_exists"]) for candidate in candidates
            ):
                reason = "invalid_canonical_target"
            elif any(bool(candidate["_is_obsolete"]) for candidate in candidates):
                reason = "obsolete_excluded"
            elif any(
                _candidate_star_rating(candidate) < self.min_star_rating
                for candidate in candidates
            ):
                reason = "below_min_star_rating"
            else:
                reason = "not_found"
            output.append(
                {
                    "GroupId": row.group_id,
                    "InputId": row.input_id,
                    "InputNamespace": self.namespace,
                    "Reason": reason,
                }
            )
        frame = pl.DataFrame(
            output,
            schema={
                "GroupId": pl.String,
                "InputId": pl.String,
                "InputNamespace": pl.String,
                "Reason": pl.Enum(sorted(_REASONS)),
            },
        )
        return self._finalize(frame)

    def extract_compounds(self) -> pl.DataFrame:
        """Return one canonical compound row per accepted match.

        Examples:
            Read the scalar annotation profile:

            >>> selection.extract_compounds().select(  # doctest: +SKIP
            ...     "ChEBIId", "PreferredName", "Formula"
            ... ).columns
            ['ChEBIId', 'PreferredName', 'Formula']
        """
        return self._extract_joined(
            """
            SELECT DISTINCT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                compound.chebi_id AS "ChEBIId",
                compound.preferred_name AS "PreferredName",
                compound.definition AS "Definition",
                compound.star_rating AS "StarRating",
                compound.is_obsolete AS "IsObsolete",
                compound.formula AS "Formula",
                compound.charge AS "Charge",
                compound.average_mass AS "AverageMass",
                compound.monoisotopic_mass AS "MonoisotopicMass",
                compound.smiles AS "SMILES",
                compound.inchi AS "InChI",
                compound.inchi_key AS "InChIKey"
            FROM _selected_compound AS selected
            JOIN compound USING (chebi_id)
            ORDER BY "GroupId", "InputId", "ChEBIId"
            """
        )

    def extract_names(self) -> pl.DataFrame:
        """Return synonyms for selected compounds.

        Examples:
            Retain each synonym's official scope:

            >>> selection.extract_names().select("Name", "Scope").columns  # doctest: +SKIP
            ['Name', 'Scope']
        """
        self._require_tables({"compound_name"}, "extract compound names")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.chebi_id AS "ChEBIId",
                names.name AS "Name",
                names.scope AS "Scope"
            FROM _selected_compound AS selected
            JOIN compound_name AS names USING (chebi_id)
            ORDER BY "GroupId", "InputId", "ChEBIId", "Scope", "Name"
            """
        )

    def extract_cross_references(self) -> pl.DataFrame:
        """Return external database accessions for selected compounds.

        Examples:
            Keep prefix and accession separate for exact reuse:

            >>> selection.extract_cross_references().select(  # doctest: +SKIP
            ...     "SourcePrefix", "Accession"
            ... ).columns
            ['SourcePrefix', 'Accession']
        """
        self._require_tables(
            {"compound_cross_reference"}, "extract compound cross-references"
        )
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.chebi_id AS "ChEBIId",
                xref.source_prefix AS "SourcePrefix",
                xref.accession AS "Accession",
                xref.xref_id AS "CrossReferenceId"
            FROM _selected_compound AS selected
            JOIN compound_cross_reference AS xref USING (chebi_id)
            ORDER BY "GroupId", "InputId", "ChEBIId", "SourcePrefix", "Accession"
            """
        )

    def extract_relations(
        self,
        *,
        direction: Literal["outgoing", "incoming", "both"] = "both",
    ) -> pl.DataFrame:
        """Return direct compound relations in the requested direction.

        Examples:
            Request only edges leaving the selected compound:

            >>> selection.extract_relations(  # doctest: +SKIP
            ...     direction="outgoing"
            ... ).columns[-3:]
            ['RelationType', 'RelationId', 'ObjectChEBIId']
        """
        if direction not in {"outgoing", "incoming", "both"}:
            raise ValueError("direction must be 'outgoing', 'incoming', or 'both'")
        self._require_tables({"compound_relation"}, "extract compound relations")
        predicates: list[str] = []
        if direction in {"outgoing", "both"}:
            predicates.append("relation.subject_chebi_id = selected.chebi_id")
        if direction in {"incoming", "both"}:
            predicates.append("relation.object_chebi_id = selected.chebi_id")
        return self._extract_joined(
            f"""
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.chebi_id AS "ChEBIId",
                relation.subject_chebi_id AS "SubjectChEBIId",
                relation.relation_type AS "RelationType",
                relation.relation_id AS "RelationId",
                relation.object_chebi_id AS "ObjectChEBIId"
            FROM _selected_compound AS selected
            JOIN compound_relation AS relation
              ON {" OR ".join(predicates)}
            ORDER BY "GroupId", "InputId", "SubjectChEBIId",
                     "RelationType", "ObjectChEBIId"
            """
        )

    def extract_structures(self) -> pl.DataFrame:
        """Return SDF molfile records for selected compounds.

        Examples:
            Keep large molfiles outside the scalar compound profile:

            >>> selection.extract_structures().select(  # doctest: +SKIP
            ...     "ChEBIId", "Molfile"
            ... ).columns
            ['ChEBIId', 'Molfile']
        """
        self._require_tables({"compound_structure"}, "extract compound structures")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.chebi_id AS "ChEBIId",
                structure.structure_index AS "StructureIndex",
                structure.molfile AS "Molfile"
            FROM _selected_compound AS selected
            JOIN compound_structure AS structure USING (chebi_id)
            ORDER BY "GroupId", "InputId", "ChEBIId", "StructureIndex"
            """
        )

    def extract_wurcs(self) -> pl.DataFrame:
        """Return WURCS representations for selected compounds.

        Examples:
            Extract source WURCS strings without parsing them:

            >>> selection.extract_wurcs().columns[-1]  # doctest: +SKIP
            'WURCS'
        """
        self._require_tables({"compound_wurcs"}, "extract compound WURCS")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "GroupId",
                selected.input_id AS "InputId",
                selected.input_namespace AS "InputNamespace",
                selected.chebi_id AS "ChEBIId",
                wurcs.wurcs AS "WURCS"
            FROM _selected_compound AS selected
            JOIN compound_wurcs AS wurcs USING (chebi_id)
            ORDER BY "GroupId", "InputId", "ChEBIId", "WURCS"
            """
        )

    def extract_ancestors(self) -> pl.DataFrame:
        """Traverse transitive ``is_a`` parents cycle-safely.

        Examples:
            Retain traversal depth for each ancestor:

            >>> selection.extract_ancestors().columns[-2:]  # doctest: +SKIP
            ['AncestorChEBIId', 'Depth']
        """
        return self._extract_hierarchy(reverse=False)

    def extract_descendants(self) -> pl.DataFrame:
        """Traverse transitive ``is_a`` children cycle-safely.

        Examples:
            Retain traversal depth for each descendant:

            >>> selection.extract_descendants().columns[-2:]  # doctest: +SKIP
            ['DescendantChEBIId', 'Depth']
        """
        return self._extract_hierarchy(reverse=True)

    def _extract_hierarchy(self, *, reverse: bool) -> pl.DataFrame:
        self._require_tables({"compound_relation"}, "traverse compound hierarchy")
        if reverse:
            first_join = "relation.object_chebi_id = selected.chebi_id"
            next_id = "relation.subject_chebi_id"
            recursive_join = "relation.object_chebi_id = hierarchy.related_id"
            label = "DescendantChEBIId"
        else:
            first_join = "relation.subject_chebi_id = selected.chebi_id"
            next_id = "relation.object_chebi_id"
            recursive_join = "relation.subject_chebi_id = hierarchy.related_id"
            label = "AncestorChEBIId"
        return self._extract_joined(
            f"""
            WITH RECURSIVE hierarchy(
                group_id, input_id, input_namespace, chebi_id,
                related_id, depth, path
            ) AS (
                SELECT
                    selected.group_id,
                    selected.input_id,
                    selected.input_namespace,
                    selected.chebi_id,
                    {next_id},
                    1,
                    [selected.chebi_id, {next_id}]
                FROM _selected_compound AS selected
                JOIN compound_relation AS relation ON {first_join}
                WHERE relation.relation_type = 'is_a'
                UNION ALL
                SELECT
                    hierarchy.group_id,
                    hierarchy.input_id,
                    hierarchy.input_namespace,
                    hierarchy.chebi_id,
                    {next_id},
                    hierarchy.depth + 1,
                    list_append(hierarchy.path, {next_id})
                FROM hierarchy
                JOIN compound_relation AS relation ON {recursive_join}
                WHERE relation.relation_type = 'is_a'
                  AND NOT list_contains(hierarchy.path, {next_id})
            )
            SELECT DISTINCT
                group_id AS "GroupId",
                input_id AS "InputId",
                input_namespace AS "InputNamespace",
                chebi_id AS "ChEBIId",
                related_id AS "{label}",
                depth AS "Depth"
            FROM hierarchy
            ORDER BY "GroupId", "InputId", "Depth", "{label}"
            """
        )

    def _candidates(self) -> pl.DataFrame:
        cached = self._candidate_cache
        if cached is not None:
            return cached.clone()
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        with duckdb.connect(str(publication.path), read_only=True) as connection:
            _create_input_table(connection, self._input_rows)
            if self.namespace == "chebi":
                source = """
                    SELECT input.group_id, input.input_id,
                           compound.chebi_id, 'primary_id' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    LEFT JOIN compound
                      ON compound.chebi_id = input.lookup_value
                    UNION ALL
                    SELECT input.group_id, input.input_id,
                           secondary.chebi_id, 'secondary_id' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    JOIN secondary_id AS secondary
                      ON secondary.secondary_chebi_id = input.lookup_value
                """
            elif self.namespace in {"inchi", "inchi_key"}:
                column = self.namespace
                source = f"""
                    SELECT input.group_id, input.input_id,
                           compound.chebi_id, '{column}' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    LEFT JOIN compound
                      ON compound.{column} = input.lookup_value
                """
            else:
                source = """
                    SELECT input.group_id, input.input_id,
                           xref.chebi_id, 'cross_reference' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    LEFT JOIN compound_cross_reference AS xref
                      ON xref.source_prefix = ?
                     AND xref.accession = input.lookup_value
                """
            parameters = (
                [] if self.namespace in _BUILTIN_NAMESPACES else [self.namespace]
            )
            rows = connection.execute(
                f"""
                WITH candidates AS ({source})
                SELECT
                    candidates.group_id,
                    candidates.input_id,
                    candidates.chebi_id,
                    candidates.match_type,
                    candidates.matched_value,
                    compound.chebi_id IS NOT NULL AS target_exists,
                    compound.is_obsolete,
                    compound.star_rating
                FROM candidates
                LEFT JOIN compound USING (chebi_id)
                """,
                parameters,
            ).fetchall()
        frame = pl.DataFrame(
            rows,
            schema={
                "GroupId": pl.String,
                "InputId": pl.String,
                "_chebi_id": pl.String,
                "_match_type": pl.String,
                "_matched_value": pl.String,
                "_target_exists": pl.Boolean,
                "_is_obsolete": pl.Boolean,
                "_star_rating": pl.Int8,
            },
            orient="row",
        ).filter(pl.col("_chebi_id").is_not_null())
        object.__setattr__(self, "_candidate_cache", frame)
        return frame.clone()

    def _invalid_target_inputs(self) -> set[tuple[str | None, str]]:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        with duckdb.connect(str(publication.path), read_only=True) as connection:
            metadata_tables = {
                row[0]
                for row in connection.execute(
                    """
                    SELECT table_name FROM information_schema.tables
                    WHERE table_schema = '_bioextract'
                    """
                ).fetchall()
            }
            if "validation_issue" not in metadata_tables:
                return set()
            _create_input_table(connection, self._input_rows)
            rows = connection.execute(
                """
                SELECT DISTINCT input.group_id, input.input_id
                FROM _input_id AS input
                JOIN _bioextract.validation_issue AS issue
                  ON issue.issue_code = 'foreign_key_violation'
                 AND issue.identifier_namespace = ?
                 AND issue.identifier_value = input.lookup_value
                """,
                [self.namespace],
            ).fetchall()
        return {(row[0], str(row[1])) for row in rows}

    def _extract_joined(self, query: str) -> pl.DataFrame:
        matches = self.extract_matches()
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        with duckdb.connect(str(publication.path), read_only=True) as connection:
            _create_selected_table(connection, matches, grouped=self._is_grouped)
            cursor = connection.execute(query)
            columns = [description[0] for description in cursor.description]
            frame = pl.DataFrame(
                cursor.fetchall(),
                schema=columns,
                orient="row",
                infer_schema_length=None,
                strict=False,
            )
        return self._finalize(frame)

    def _require_tables(self, required: set[str], operation: str) -> None:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        missing = sorted(required - set(publication.tables))
        if missing:
            raise ChEBICapabilityError(
                f"ChEBI publication cannot {operation}; missing relations: {missing}"
            )

    def _finalize(self, frame: pl.DataFrame) -> pl.DataFrame:
        if not self._is_grouped and "GroupId" in frame.columns:
            return frame.drop("GroupId")
        return frame


def open_chebi_publication(path: str | Path) -> _ChEBIPublication:
    """Validate and describe one bioextract ChEBI DuckDB publication."""
    publication_path = Path(path)
    if not publication_path.exists():
        raise FileNotFoundError(publication_path)
    if not publication_path.is_file():
        raise ValueError(f"ChEBI DuckDB path is not a file: {publication_path}")
    try:
        connection = duckdb.connect(str(publication_path), read_only=True)
    except duckdb.Error as error:
        raise ValueError(
            f"Cannot open ChEBI DuckDB publication: {publication_path}"
        ) from error
    try:
        metadata_tables = {
            str(row[0])
            for row in connection.execute(
                """
                SELECT table_name FROM information_schema.tables
                WHERE table_schema = '_bioextract'
                  AND table_type = 'BASE TABLE'
                """
            ).fetchall()
        }
        missing = sorted(_BASE_METADATA_TABLES - metadata_tables)
        if missing:
            raise ValueError(
                f"DuckDB file is missing bioextract metadata tables: {missing}"
            )
        metadata = {
            str(key): str(value)
            for key, value in connection.execute(
                "SELECT key, value FROM _bioextract.metadata"
            ).fetchall()
        }
        metadata_version = metadata.get("bioextract.metadata_schema_version")
        if metadata_version not in _METADATA_VERSIONS:
            raise ValueError(
                f"Unsupported ChEBI metadata schema version: {metadata_version!r}"
            )
        if metadata_version == "3":
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
                raise ValueError(f"ChEBI metadata v3 is missing keys: {missing_v3}")
            if (
                metadata.get("bioextract.source_schema_profile")
                != SOURCE_SCHEMA_PROFILE
            ):
                raise ValueError("Unsupported ChEBI source schema profile")
        if metadata_version == "2":
            validate_duckdb_validation_state(connection, metadata)
        if metadata.get("bioextract.resource_name") != "chebi":
            raise ValueError("DuckDB file is not a bioextract ChEBI publication")
        resource_schema_version = metadata.get(
            "bioextract.resource_schema_version"
            if metadata_version == "3"
            else "bioextract.schema_version"
        )
        if resource_schema_version != "chebi-duckdb-v1":
            raise ValueError("Unsupported ChEBI resource schema version")
        tables = frozenset(
            str(row[0])
            for row in connection.execute(
                """
                SELECT table_name FROM information_schema.tables
                WHERE table_schema = 'main' AND table_type = 'BASE TABLE'
                """
            ).fetchall()
        )
        required = {
            "compound",
            "secondary_id",
            "compound_cross_reference",
            "compound_relation",
        }
        missing = sorted(required - set(tables))
        if missing:
            raise ValueError(f"ChEBI publication lacks canonical relations: {missing}")
        recorded = {
            str(name): int(count)
            for name, count in connection.execute(
                "SELECT table_name, row_count FROM _bioextract.table_info"
            ).fetchall()
        }
        if set(recorded) != set(tables):
            raise ValueError("ChEBI publication table inventory does not match")
        for table_name, expected in recorded.items():
            observed_row = connection.execute(
                f'SELECT count(*) FROM "{table_name}"'
            ).fetchone()
            if observed_row is None:
                raise ValueError(f"Cannot count ChEBI publication table {table_name}")
            observed = int(observed_row[0])
            if observed != expected:
                raise ValueError(
                    f"ChEBI publication row count mismatch for {table_name}"
                )
        chebi_type = connection.execute(
            """
            SELECT data_type FROM information_schema.columns
            WHERE table_schema='main' AND table_name='compound'
              AND column_name='chebi_id'
            """
        ).fetchone()
        if chebi_type is None or str(chebi_type[0]) != "VARCHAR":
            raise ValueError("compound.chebi_id must be a VARCHAR CURIE")
        malformed_row = connection.execute(
            """
            SELECT count(*) FROM compound
            WHERE NOT regexp_full_match(chebi_id, 'CHEBI:[0-9]+')
            """
        ).fetchone()
        if malformed_row is None:
            raise ValueError("Cannot validate ChEBI compound identifiers")
        malformed = int(malformed_row[0])
        if malformed:
            raise ValueError("compound.chebi_id contains malformed CURIE values")
        namespaces = set(_BUILTIN_NAMESPACES)
        namespaces.update(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT source_prefix FROM compound_cross_reference"
            ).fetchall()
            if row[0]
        )
    except duckdb.Error as error:
        raise ValueError(
            f"Invalid bioextract ChEBI publication: {publication_path}"
        ) from error
    finally:
        connection.close()
    return _ChEBIPublication(
        path=publication_path,
        tables=tables,
        metadata=metadata,
        namespaces=frozenset(namespaces),
    )


def create_selection(
    database: ChEBIDatabase,
    ids: Iterable[str],
    *,
    namespace: str,
    min_star_rating: int,
    include_obsolete: bool,
) -> ChEBICompoundSelection:
    return _selection(
        database,
        ((None, value) for value in ids),
        namespace=namespace,
        min_star_rating=min_star_rating,
        include_obsolete=include_obsolete,
        grouped=False,
    )


def create_group_selection(
    database: ChEBIDatabase,
    ids_by_group: Mapping[str, Iterable[str]],
    *,
    namespace: str,
    min_star_rating: int,
    include_obsolete: bool,
) -> ChEBICompoundSelection:
    return _selection(
        database,
        (
            (str(group_id), value)
            for group_id, values in ids_by_group.items()
            for value in values
        ),
        namespace=namespace,
        min_star_rating=min_star_rating,
        include_obsolete=include_obsolete,
        grouped=True,
    )


def _selection(
    database: ChEBIDatabase,
    inputs: Iterable[tuple[str | None, str]],
    *,
    namespace: str,
    min_star_rating: int,
    include_obsolete: bool,
    grouped: bool,
) -> ChEBICompoundSelection:
    publication = database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
    normalized_namespace = str(namespace).strip().lower()
    if normalized_namespace not in publication.namespaces:
        available = ", ".join(sorted(publication.namespaces))
        raise ValueError(
            f"Unknown ChEBI namespace: {namespace!r}; available: {available}"
        )
    if not 1 <= min_star_rating <= 3:
        raise ValueError("min_star_rating must be between 1 and 3")
    rows: set[tuple[str | None, str, str]] = set()
    for group_id, value in inputs:
        text = str(value).strip()
        if not text:
            continue
        lookup = _normalize_lookup(text, normalized_namespace)
        rows.add(
            (group_id, lookup if normalized_namespace == "chebi" else text, lookup)
        )
    input_rows = tuple(
        _InputRow(group_id, input_id, lookup)
        for group_id, input_id, lookup in sorted(
            rows, key=lambda row: (row[0] or "", row[1], row[2])
        )
    )
    return ChEBICompoundSelection(
        database=database,
        namespace=normalized_namespace,
        min_star_rating=min_star_rating,
        include_obsolete=bool(include_obsolete),
        _input_rows=input_rows,
        _is_grouped=grouped,
    )


def _normalize_lookup(value: str, namespace: str) -> str:
    if namespace != "chebi":
        return value
    match = _CHEBI_ID.fullmatch(value)
    if match is None:
        raise ValueError(f"Invalid ChEBI identifier: {value!r}")
    return f"CHEBI:{int(match.group(1))}"


def _candidate_star_rating(candidate: Mapping[str, object]) -> int:
    value = candidate.get("_star_rating")
    return 0 if value is None else int(str(value))


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
        CREATE TEMP TABLE _selected_compound (
            group_id VARCHAR,
            input_id VARCHAR NOT NULL,
            input_namespace VARCHAR NOT NULL,
            chebi_id VARCHAR NOT NULL
        )
        """
    )
    rows = [
        (
            row["GroupId"] if grouped else None,
            row["InputId"],
            row["InputNamespace"],
            row["ChEBIId"],
        )
        for row in matches.to_dicts()
    ]
    if rows:
        connection.executemany(
            "INSERT INTO _selected_compound VALUES (?, ?, ?, ?)",
            rows,
        )
