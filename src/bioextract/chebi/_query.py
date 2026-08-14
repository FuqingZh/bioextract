from __future__ import annotations

import copy
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import register_deferred_frame_source
from bioextract._publication import validate_duckdb_metadata_v2
from bioextract._shared import validate_group_ids
from bioextract.errors import CapabilityError

if TYPE_CHECKING:
    from .chebi import ChEBIDatabase

_METADATA_VERSIONS = {"2"}
SOURCE_SCHEMA_PROFILE = "chebi-release-bundle-v1"
_BASE_METADATA_TABLES = {
    "metadata",
    "source_file",
    "table_info",
    "column_mapping",
    "validation_issue",
}
_CHEBI_ID = re.compile(r"^(?:CHEBI:)?([0-9]+)$", re.IGNORECASE)
_BUILTIN_NAMESPACES = {"chebi", "inchi", "inchi_key"}
_REASONS = {
    "not_found",
    "below_min_star_rating",
    "obsolete_excluded",
    "invalid_canonical_target",
}
_CHEBI_RELATION_TYPES: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
    "input_namespace": pl.String,
    "chebi_id": pl.String,
    "match_type": pl.String,
    "matched_value": pl.String,
    "reason": pl.Enum(sorted(_REASONS)),
    "preferred_name": pl.String,
    "definition": pl.String,
    "star_rating": pl.Int8,
    "is_obsolete": pl.Boolean,
    "formula": pl.String,
    "charge": pl.Int32,
    "average_mass": pl.Float64,
    "monoisotopic_mass": pl.Float64,
    "smiles": pl.String,
    "inchi": pl.String,
    "inchi_key": pl.String,
    "name": pl.String,
    "scope": pl.String,
    "source_prefix": pl.String,
    "accession": pl.String,
    "xref_id": pl.String,
    "subject_chebi_id": pl.String,
    "relation_type": pl.String,
    "relation_id": pl.String,
    "object_chebi_id": pl.String,
    "structure_index": pl.Int32,
    "molfile": pl.String,
    "wurcs": pl.String,
    "ancestor_chebi_id": pl.String,
    "descendant_chebi_id": pl.String,
    "depth": pl.Int64,
}


def _relation_schema(columns: Iterable[str]) -> SchemaDict:
    return {column: _CHEBI_RELATION_TYPES[column] for column in columns}


@dataclass(frozen=True, slots=True)
class _ChEBIPublication:
    path: Path
    tables: frozenset[str]
    metadata: Mapping[str, str]
    namespaces: frozenset[str]


@dataclass(frozen=True, slots=True)
class _InputRow:
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
    _group_ids: tuple[str, ...] = field(repr=False)
    _group_membership: tuple[tuple[str, str], ...] = field(repr=False)
    _is_grouped: bool = field(repr=False)
    _candidate_cache: pl.DataFrame | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _invalid_target_cache: frozenset[str] | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def matches(self) -> pl.LazyFrame:
        """Return identifier-to-canonical-compound matches lazily.

        Examples:
            >>> selection.matches().select("chebi_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "match_type",
            "matched_value",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_matches,
        )

    def _eager_matches(self) -> pl.DataFrame:
        """Return identifier-to-canonical-compound matches.

        Examples:
            Observe the canonical shared identifier:

            >>> selection.matches().select("chebi_id").collect()  # doctest: +SKIP
            ['CHEBI:15377']
        """
        candidates = self._candidates()
        if candidates.is_empty():
            frame = pl.DataFrame(
                schema={
                    "input_id": pl.String,
                    "input_namespace": pl.String,
                    "chebi_id": pl.String,
                    "match_type": pl.String,
                    "matched_value": pl.String,
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
                    "input_id",
                    pl.lit(self.namespace).alias("input_namespace"),
                    pl.col("_chebi_id").alias("chebi_id"),
                    pl.col("_match_type").alias("match_type"),
                    pl.col("_matched_value").alias("matched_value"),
                )
                .unique()
                .sort(["input_id", "chebi_id"])
            )
        frame = self._expand_groups(frame)
        sort_cols = (
            ["group_id", "input_id", "chebi_id"]
            if self._is_grouped
            else [
                "input_id",
                "chebi_id",
            ]
        )
        return frame.sort(sort_cols)

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return policy-rejected inputs and stable reasons lazily.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 3)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "reason",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_unmatched_ids,
        )

    def _eager_unmatched_ids(self) -> pl.DataFrame:
        """Return inputs without a policy-accepted match and a stable reason.

        Examples:
            Inspect stable reason codes separately from successful matches:

            >>> unmatched = selection.unmatched_ids().collect()  # doctest: +SKIP
            >>> "reason" in unmatched.columns  # doctest: +SKIP
            True
        """
        matches = set(self._eager_matches().get_column("input_id").to_list())
        candidate_rows = self._candidates().to_dicts()
        by_input: dict[str, list[dict[str, object]]] = {}
        for row in candidate_rows:
            by_input.setdefault(str(row["input_id"]), []).append(row)
        issue_inputs = self._invalid_target_inputs()
        output: list[dict[str, str]] = []
        for row in self._input_rows:
            if row.input_id in matches:
                continue
            candidates = by_input.get(row.input_id, [])
            if row.input_id in issue_inputs or any(
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
                    "input_id": row.input_id,
                    "input_namespace": self.namespace,
                    "reason": reason,
                }
            )
        frame = pl.DataFrame(
            output,
            schema={
                "input_id": pl.String,
                "input_namespace": pl.String,
                "reason": pl.Enum(sorted(_REASONS)),
            },
        )
        frame = self._expand_groups(frame)
        sort_cols = ["group_id", "input_id"] if self._is_grouped else ["input_id"]
        return frame.sort(sort_cols)

    def compounds(self) -> pl.LazyFrame:
        """Return one canonical compound row per accepted match lazily.

        Examples:
            >>> selection.compounds().select("chebi_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "preferred_name",
            "definition",
            "star_rating",
            "is_obsolete",
            "formula",
            "charge",
            "average_mass",
            "monoisotopic_mass",
            "smiles",
            "inchi",
            "inchi_key",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_compounds,
        )

    def _eager_compounds(self) -> pl.DataFrame:
        """Return one canonical compound row per accepted match.

        Examples:
            Read the scalar annotation profile:

            >>> selection.compounds().select(  # doctest: +SKIP
            ...     "chebi_id", "preferred_name", "formula"
            ... ).collect_schema().names()
            ['chebi_id', 'preferred_name', 'formula']
        """
        return self._extract_joined(
            """
            SELECT DISTINCT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                compound.chebi_id AS chebi_id,
                compound.preferred_name AS preferred_name,
                compound.definition AS definition,
                compound.star_rating AS star_rating,
                compound.is_obsolete AS is_obsolete,
                compound.formula AS formula,
                compound.charge AS charge,
                compound.average_mass AS average_mass,
                compound.monoisotopic_mass AS monoisotopic_mass,
                compound.smiles AS smiles,
                compound.inchi AS inchi,
                compound.inchi_key AS inchi_key
            FROM _selected_compound AS selected
            JOIN compound USING (chebi_id)
            ORDER BY group_id, input_id, chebi_id
            """
        )

    def names(self) -> pl.LazyFrame:
        """Return synonyms for selected compounds lazily.

        Examples:
            >>> selection.names().select("name", "scope").collect()  # doctest: +SKIP
            shape: (..., 2)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "name",
            "scope",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_names,
        )

    def _eager_names(self) -> pl.DataFrame:
        """Return synonyms for selected compounds.

        Examples:
            Retain each synonym's official scope:

            >>> selection.names().select("name", "scope").collect_schema().names()  # doctest: +SKIP
            ['name', 'scope']
        """
        self._require_tables({"compound_name"}, "extract compound names")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.chebi_id AS chebi_id,
                names.name AS name,
                names.scope AS scope
            FROM _selected_compound AS selected
            JOIN compound_name AS names USING (chebi_id)
            ORDER BY group_id, input_id, chebi_id, scope, name
            """
        )

    def cross_references(self) -> pl.LazyFrame:
        """Return external database accessions for selected compounds lazily.

        Examples:
            >>> selection.cross_references().collect()  # doctest: +SKIP
            shape: (..., 7)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "source_prefix",
            "accession",
            "xref_id",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_cross_references,
        )

    def _eager_cross_references(self) -> pl.DataFrame:
        """Return external database accessions for selected compounds.

        Examples:
            Keep prefix and accession separate for exact reuse:

            >>> selection.cross_references().select(  # doctest: +SKIP
            ...     "source_prefix", "accession"
            ... ).collect_schema().names()
            ['source_prefix', 'accession']
        """
        self._require_tables(
            {"compound_cross_reference"}, "extract compound cross-references"
        )
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.chebi_id AS chebi_id,
                xref.source_prefix AS source_prefix,
                xref.accession AS accession,
                xref.xref_id AS xref_id
            FROM _selected_compound AS selected
            JOIN compound_cross_reference AS xref USING (chebi_id)
            ORDER BY group_id, input_id, chebi_id, source_prefix, accession
            """
        )

    def relations(
        self,
        *,
        direction: Literal["outgoing", "incoming", "both"] = "both",
    ) -> pl.LazyFrame:
        """Return direct compound relations in the requested direction lazily.

        Examples:
            >>> selection.relations(direction="outgoing").collect()  # doctest: +SKIP
            shape: (..., 8)
        """
        if direction not in {"outgoing", "incoming", "both"}:
            raise ValueError("direction must be 'outgoing', 'incoming', or 'both'")
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "subject_chebi_id",
            "relation_type",
            "relation_id",
            "object_chebi_id",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=lambda: snapshot._eager_relations(direction=direction),
        )

    def _eager_relations(
        self,
        *,
        direction: Literal["outgoing", "incoming", "both"] = "both",
    ) -> pl.DataFrame:
        """Return direct compound relations in the requested direction.

        Examples:
            Request only edges leaving the selected compound:

            >>> selection.relations(  # doctest: +SKIP
            ...     direction="outgoing"
            ... ).collect_schema().names()[-3:]
            ['relation_type', 'relation_id', 'object_chebi_id']
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
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.chebi_id AS chebi_id,
                relation.subject_chebi_id AS subject_chebi_id,
                relation.relation_type AS relation_type,
                relation.relation_id AS relation_id,
                relation.object_chebi_id AS object_chebi_id
            FROM _selected_compound AS selected
            JOIN compound_relation AS relation
              ON {" OR ".join(predicates)}
            ORDER BY group_id, input_id, subject_chebi_id,
                     relation_type, object_chebi_id
            """
        )

    def structures(self) -> pl.LazyFrame:
        """Return SDF molfile records for selected compounds lazily.

        Examples:
            >>> selection.structures().select("chebi_id", "molfile").collect()  # doctest: +SKIP
            shape: (..., 2)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "structure_index",
            "molfile",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_structures,
        )

    def _eager_structures(self) -> pl.DataFrame:
        """Return SDF molfile records for selected compounds.

        Examples:
            Keep large molfiles outside the scalar compound profile:

            >>> selection.structures().select(  # doctest: +SKIP
            ...     "chebi_id", "molfile"
            ... ).collect_schema().names()
            ['chebi_id', 'molfile']
        """
        self._require_tables({"compound_structure"}, "extract compound structures")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.chebi_id AS chebi_id,
                structure.structure_index AS structure_index,
                structure.molfile AS molfile
            FROM _selected_compound AS selected
            JOIN compound_structure AS structure USING (chebi_id)
            ORDER BY group_id, input_id, chebi_id, structure_index
            """
        )

    def wurcs(self) -> pl.LazyFrame:
        """Return WURCS representations for selected compounds lazily.

        Examples:
            >>> selection.wurcs().select("chebi_id", "wurcs").collect()  # doctest: +SKIP
            shape: (..., 2)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "wurcs",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=snapshot._eager_wurcs,
        )

    def _eager_wurcs(self) -> pl.DataFrame:
        """Return WURCS representations for selected compounds.

        Examples:
            Extract source WURCS strings without parsing them:

            >>> selection.wurcs().collect_schema().names()[-1]  # doctest: +SKIP
            'wurcs'
        """
        self._require_tables({"compound_wurcs"}, "extract compound WURCS")
        return self._extract_joined(
            """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.chebi_id AS chebi_id,
                wurcs.wurcs AS wurcs
            FROM _selected_compound AS selected
            JOIN compound_wurcs AS wurcs USING (chebi_id)
            ORDER BY group_id, input_id, chebi_id, wurcs
            """
        )

    def ancestors(self) -> pl.LazyFrame:
        """Traverse transitive ``is_a`` parents cycle-safely and lazily.

        Examples:
            >>> selection.ancestors().select("ancestor_chebi_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "ancestor_chebi_id",
            "depth",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=lambda: snapshot._extract_hierarchy(reverse=False),
        )

    def descendants(self) -> pl.LazyFrame:
        """Traverse transitive ``is_a`` children cycle-safely and lazily.

        Examples:
            >>> selection.descendants().select("descendant_chebi_id").collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        snapshot = copy.copy(self)
        columns = (["group_id"] if self._is_grouped else []) + [
            "input_id",
            "input_namespace",
            "chebi_id",
            "descendant_chebi_id",
            "depth",
        ]
        return register_deferred_frame_source(
            schema=_relation_schema(columns),
            frame=lambda: snapshot._extract_hierarchy(reverse=True),
        )

    def _extract_hierarchy(self, *, reverse: bool) -> pl.DataFrame:
        self._require_tables({"compound_relation"}, "traverse compound hierarchy")
        if reverse:
            first_join = "relation.object_chebi_id = selected.chebi_id"
            next_id = "relation.subject_chebi_id"
            recursive_join = "relation.object_chebi_id = hierarchy.related_id"
            label = "descendant_chebi_id"
        else:
            first_join = "relation.subject_chebi_id = selected.chebi_id"
            next_id = "relation.object_chebi_id"
            recursive_join = "relation.subject_chebi_id = hierarchy.related_id"
            label = "ancestor_chebi_id"
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
                group_id AS "group_id",
                input_id AS "input_id",
                input_namespace AS "input_namespace",
                chebi_id AS chebi_id,
                related_id AS "{label}",
                depth AS depth
            FROM hierarchy
            ORDER BY group_id, input_id, depth, "{label}"
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
                    SELECT input.input_id,
                           compound.chebi_id, 'primary_id' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    LEFT JOIN compound
                      ON compound.chebi_id = input.lookup_value
                    UNION ALL
                    SELECT input.input_id,
                           secondary.chebi_id, 'secondary_id' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    JOIN secondary_id AS secondary
                      ON secondary.secondary_chebi_id = input.lookup_value
                """
            elif self.namespace in {"inchi", "inchi_key"}:
                column = self.namespace
                source = f"""
                    SELECT input.input_id,
                           compound.chebi_id, '{column}' AS match_type,
                           input.lookup_value AS matched_value
                    FROM _input_id AS input
                    LEFT JOIN compound
                      ON compound.{column} = input.lookup_value
                """
            else:
                source = """
                    SELECT input.input_id,
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
                "input_id": pl.String,
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

    def _invalid_target_inputs(self) -> frozenset[str]:
        cached = self._invalid_target_cache
        if cached is not None:
            return cached
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
                result = frozenset[str]()
                object.__setattr__(self, "_invalid_target_cache", result)
                return result
            _create_input_table(connection, self._input_rows)
            rows = connection.execute(
                """
                SELECT DISTINCT input.input_id
                FROM _input_id AS input
                JOIN _bioextract.validation_issue AS issue
                  ON issue.issue_code = 'foreign_key_violation'
                 AND issue.identifier_namespace = ?
                 AND issue.identifier_value = input.lookup_value
                """,
                [self.namespace],
            ).fetchall()
        result = frozenset(str(row[0]) for row in rows)
        object.__setattr__(self, "_invalid_target_cache", result)
        return result

    def _extract_joined(self, query: str) -> pl.DataFrame:
        matches = self._eager_matches()
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
            raise CapabilityError(
                f"ChEBI publication cannot {operation}; missing relations: {missing}"
            )

    def _expand_groups(self, frame: pl.DataFrame) -> pl.DataFrame:
        if not self._is_grouped:
            return frame
        membership = pl.DataFrame(
            self._group_membership,
            schema={"group_id": pl.String, "input_id": pl.String},
            orient="row",
        )
        return membership.join(frame, on="input_id", how="inner").select(
            "group_id", *frame.columns
        )

    def _finalize(self, frame: pl.DataFrame) -> pl.DataFrame:
        if not self._is_grouped and "group_id" in frame.columns:
            return frame.drop("group_id")
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
        validate_duckdb_metadata_v2(connection, metadata)
        if metadata.get("bioextract.source_schema_profile") != SOURCE_SCHEMA_PROFILE:
            raise ValueError("Unsupported ChEBI source schema profile")
        if metadata.get("bioextract.resource_name") != "chebi":
            raise ValueError("DuckDB file is not a bioextract ChEBI publication")
        resource_schema_version = metadata.get("bioextract.resource_schema_version")
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
        group_ids=(),
    )


def create_group_selection(
    database: ChEBIDatabase,
    ids_by_group: Mapping[str, Iterable[str]],
    *,
    namespace: str,
    min_star_rating: int,
    include_obsolete: bool,
) -> ChEBICompoundSelection:
    group_ids: list[str] = []
    inputs: list[tuple[str, str]] = []
    for raw_group_id, values in ids_by_group.items():
        group_id = str(raw_group_id).strip()
        group_ids.append(group_id)
        inputs.extend((group_id, value) for value in values)
    validate_group_ids(group_ids)
    return _selection(
        database,
        inputs,
        namespace=namespace,
        min_star_rating=min_star_rating,
        include_obsolete=include_obsolete,
        grouped=True,
        group_ids=tuple(sorted(group_ids)),
    )


def _selection(
    database: ChEBIDatabase,
    inputs: Iterable[tuple[str | None, str]],
    *,
    namespace: str,
    min_star_rating: int,
    include_obsolete: bool,
    grouped: bool,
    group_ids: tuple[str, ...],
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
    rows: set[tuple[str, str]] = set()
    group_membership: set[tuple[str, str]] = set()
    for group_id, value in inputs:
        text = str(value).strip()
        if not text:
            continue
        lookup = _normalize_lookup(text, normalized_namespace)
        input_id = lookup if normalized_namespace == "chebi" else text
        rows.add((input_id, lookup))
        if grouped:
            if group_id is None:
                raise AssertionError("grouped ChEBI input lacks group_id")
            group_membership.add((group_id, input_id))
    input_rows = tuple(
        _InputRow(input_id, lookup)
        for input_id, lookup in sorted(rows, key=lambda row: (row[0], row[1]))
    )
    return ChEBICompoundSelection(
        database=database,
        namespace=normalized_namespace,
        min_star_rating=min_star_rating,
        include_obsolete=bool(include_obsolete),
        _input_rows=input_rows,
        _group_ids=group_ids,
        _group_membership=tuple(sorted(group_membership)),
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
            input_id VARCHAR NOT NULL,
            lookup_value VARCHAR NOT NULL
        )
        """
    )
    if rows:
        connection.executemany(
            "INSERT INTO _input_id VALUES (?, ?)",
            [(row.input_id, row.lookup_value) for row in rows],
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
            row["group_id"] if grouped else None,
            row["input_id"],
            row["input_namespace"],
            row["chebi_id"],
        )
        for row in matches.to_dicts()
    ]
    if rows:
        connection.executemany(
            "INSERT INTO _selected_compound VALUES (?, ?, ?, ?)",
            rows,
        )
