from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import register_replayable_source
from bioextract._publication import validate_duckdb_metadata_v1
from bioextract._shared import validate_group_ids
from bioextract.errors import CapabilityError, IntegrityError

from .constant import SCHEMA_VERSION, SOURCE_SCHEMA_PROFILE, RheaNamespace

if TYPE_CHECKING:
    from .rhea import RheaDatabase


_METADATA_TABLES = {
    "metadata",
    "source_file",
    "table_info",
    "column_mapping",
}
_METADATA_SCHEMA_VERSIONS = {"1"}
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
    "group_id": pl.String,
    "input_id": pl.String,
    "input_namespace": pl.String,
    "rhea_id": pl.Int64,
    "master_id": pl.Int64,
    "direction": pl.String,
}
_SCHEMA_UNIQUE_MATCH: SchemaDict = {
    name: dtype for name, dtype in _SCHEMA_MATCH.items() if name != "group_id"
}
_SCHEMA_REACTION: SchemaDict = {
    **_SCHEMA_MATCH,
    "accession": pl.String,
    "equation": pl.String,
    "equation_html": pl.String,
    "status": pl.String,
    "is_balanced": pl.Boolean,
    "is_transport": pl.Boolean,
    "comment": pl.String,
    "is_obsolete": pl.Boolean,
    "reaction_smiles": pl.String,
}
_SCHEMA_PARTICIPANT: SchemaDict = {
    **_SCHEMA_MATCH,
    "participant_id": pl.String,
    "compound_id": pl.String,
    "side": pl.String,
    "directional_role": pl.String,
    "coefficient_text": pl.String,
    "coefficient_numeric": pl.Float64,
    "location": pl.String,
    "rhea_compound_id": pl.Int64,
    "compound_accession": pl.String,
    "compound_type": pl.String,
    "compound_name": pl.String,
    "formula": pl.String,
    "charge_text": pl.String,
    "charge_numeric": pl.Int64,
    "chebi_id": pl.String,
    "underlying_chebi_id": pl.String,
    "polymerization_index": pl.String,
}
_SCHEMA_CROSS_REFERENCE: SchemaDict = {
    **_SCHEMA_MATCH,
    "reference_database": pl.String,
    "reference_id": pl.String,
}
_SCHEMA_UNIPROT_MAPPING: SchemaDict = {
    "rhea_id": pl.Int64,
    "master_id": pl.Int64,
    "direction": pl.String,
    "uniprot_id": pl.String,
    "uniprot_section": pl.String,
}
_SCHEMA_PUBLICATION: SchemaDict = {
    **_SCHEMA_MATCH,
    "pubmed_id": pl.String,
}
_SCHEMA_RELATIONSHIP: SchemaDict = {
    **_SCHEMA_MATCH,
    "from_reaction_id": pl.Int64,
    "to_reaction_id": pl.Int64,
    "relation_type": pl.String,
}
_SCHEMA_UNMATCHED: SchemaDict = {
    "group_id": pl.String,
    "input_id": pl.String,
    "input_namespace": pl.String,
}
_SCHEMA_UNIQUE_UNMATCHED: SchemaDict = {
    name: dtype for name, dtype in _SCHEMA_UNMATCHED.items() if name != "group_id"
}
_NESTED_OUTPUT_BATCH_SIZE = 256


def _nested_neighborhood_schema(*, grouped: bool) -> SchemaDict:
    input_fields: SchemaDict = {}
    if grouped:
        input_fields["group_id"] = pl.String
    input_fields.update(
        {
            "input_id": pl.String,
            "input_namespace": pl.String,
        }
    )
    return {
        "rhea_id": pl.Int64,
        "master_id": pl.Int64,
        "direction": pl.String,
        "inputs": pl.List(pl.Struct(input_fields)),
        "uniprot_entries": pl.List(
            pl.Struct(
                {
                    "uniprot_id": pl.String,
                    "uniprot_section": pl.String,
                }
            )
        ),
    }


@dataclass(frozen=True, slots=True)
class _RheaPublication:
    path: Path
    identity: tuple[int, int, int, int, int]
    tables: frozenset[str]
    views: frozenset[str]
    metadata: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class _InputRow:
    input_id: str
    lookup_value: str


@dataclass(frozen=True, slots=True)
class _RheaLazyPlan:
    """Immutable values required to replay one Rhea relation execution."""

    path: Path
    identity: tuple[int, int, int, int, int]
    tables: frozenset[str]
    views: frozenset[str]
    namespace: RheaNamespace
    include_obsolete: bool
    input_rows: tuple[_InputRow, ...]
    group_membership: tuple[tuple[str, str], ...]
    grouped: bool


@dataclass(frozen=True, slots=True)
class RheaReactionSelection:
    """Deferred Rhea reaction selection.

    The noun-based relation methods return native ``polars.LazyFrame``
    objects. Each frame captures immutable selection values and opens a fresh
    read-only DuckDB connection per Polars execution; eager work is private.

    Examples:
        Inspect the declared namespace without executing the query:

        >>> selection.namespace  # doctest: +SKIP
        'chebi'
    """

    database: RheaDatabase = field(repr=False)
    namespace: RheaNamespace
    include_obsolete: bool
    _input_rows: tuple[_InputRow, ...] = field(repr=False)
    _group_ids: tuple[str, ...] = field(repr=False)
    _group_membership: tuple[tuple[str, str], ...] = field(
        repr=False,
    )
    _is_grouped: bool = field(repr=False)
    _matches_cache: pl.DataFrame | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )

    def _eager_matches(self) -> pl.DataFrame:
        """Resolve input identifiers to exact Rhea reactions.

        Examples:
            Materialize the identifier-to-reaction mapping:

            >>> selection.matches().select("rhea_id").head(1).collect().to_dicts()  # doctest: +SKIP
            [10000]
        """
        frame = self._matches_cache
        if frame is None:
            frame = self._query_matches()
            object.__setattr__(self, "_matches_cache", frame)
        return frame.clone()

    def matches(self) -> pl.LazyFrame:
        """Return the normalized input-to-reaction relation lazily.

        The returned frame is replayable and does not retain an open DuckDB
        connection. Callers choose native Polars execution:

        Examples:
            >>> lf = selection.matches()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        plan = self._lazy_plan()
        return _relation_frame(
            plan,
            query=_matches_query(plan),
            raw_schema=_SCHEMA_MATCH,
            output_schema=_public_schema(_SCHEMA_MATCH, grouped=plan.grouped),
            prepare_selected=False,
        )

    def reactions(self) -> pl.LazyFrame:
        """Return selected reaction facts as a native lazy relation.

        Use ``collect()``, ``collect_batches()``, or any native Polars
        ``sink_*`` method on the returned frame.

        Examples:
            >>> lf = selection.reactions()  # doctest: +SKIP
            >>> lf.select("rhea_id", "direction").collect()  # doctest: +SKIP
            shape: (..., 2)
        """

        plan = self._lazy_plan()
        _require_tables(plan, {"reaction"}, operation="read reactions")
        has_smiles = "reaction_smiles" in plan.tables
        smiles_expression = (
            "smiles.reaction_smiles" if has_smiles else "CAST(NULL AS VARCHAR)"
        )
        smiles_join = (
            "LEFT JOIN reaction_smiles AS smiles USING (rhea_id)" if has_smiles else ""
        )
        query = f"""
            SELECT
                selected.group_id AS group_id,
                selected.input_id AS input_id,
                selected.input_namespace AS input_namespace,
                reaction.rhea_id AS rhea_id,
                reaction.master_id AS master_id,
                reaction.direction AS direction,
                reaction.accession AS accession,
                reaction.equation AS equation,
                reaction.equation_html AS equation_html,
                reaction.status AS status,
                reaction.is_balanced AS is_balanced,
                reaction.is_transport AS is_transport,
                reaction.comment AS comment,
                reaction.is_obsolete AS is_obsolete,
                {smiles_expression} AS reaction_smiles
            FROM _selected_reaction AS selected
            JOIN reaction USING (rhea_id)
            {smiles_join}
            ORDER BY group_id NULLS FIRST, input_id, rhea_id
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_REACTION,
            output_schema=_public_schema(_SCHEMA_REACTION, grouped=plan.grouped),
        )

    def participants(self) -> pl.LazyFrame:
        """Return direction-aware reaction participants lazily.

        Examples:
            >>> lf = selection.participants()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        plan = self._lazy_plan()
        _require_tables(
            plan,
            {"reaction_participant", "compound"},
            operation="read reaction participants",
        )
        if "reaction_participant_direction" not in plan.views:
            raise CapabilityError(
                "Rhea publication lacks view required to read reaction "
                "participants: reaction_participant_direction"
            )
        query = """
            SELECT
                selected.group_id AS group_id,
                selected.input_id AS input_id,
                selected.input_namespace AS input_namespace,
                participant.rhea_id AS rhea_id,
                participant.master_id AS master_id,
                participant.direction AS direction,
                participant.participant_id AS participant_id,
                participant.compound_id AS compound_id,
                participant.side AS side,
                participant.directional_role AS directional_role,
                participant.coefficient_text AS coefficient_text,
                participant.coefficient_numeric AS coefficient_numeric,
                participant.location AS location,
                compound.rhea_compound_id AS rhea_compound_id,
                compound.public_accession AS compound_accession,
                compound.compound_type AS compound_type,
                compound.name AS compound_name,
                compound.formula AS formula,
                compound.charge_text AS charge_text,
                compound.charge_numeric AS charge_numeric,
                compound.chebi_id AS chebi_id,
                compound.underlying_chebi_id AS underlying_chebi_id,
                compound.polymerization_index AS polymerization_index
            FROM _selected_reaction AS selected
            JOIN reaction_participant_direction AS participant USING (rhea_id)
            LEFT JOIN compound USING (compound_id)
            ORDER BY group_id NULLS FIRST, input_id, rhea_id, side, participant_id
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_PARTICIPANT,
            output_schema=_public_schema(_SCHEMA_PARTICIPANT, grouped=plan.grouped),
        )

    def cross_references(self) -> pl.LazyFrame:
        """Return external reaction xrefs without implicit UniProt unioning.

        UniProt rows are exposed independently by :meth:`uniprot_mappings`.
        This keeps external database identifiers and protein mappings from
        acquiring a misleading common row shape.

        Examples:
            >>> lf = selection.cross_references()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        plan = self._lazy_plan()
        _require_tables(plan, {"reaction_xref"}, operation="read cross-references")
        query = """
            SELECT DISTINCT
                selected.group_id AS group_id,
                selected.input_id AS input_id,
                selected.input_namespace AS input_namespace,
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                xref.database_name AS reference_database,
                xref.external_id AS reference_id
            FROM _selected_reaction AS selected
            JOIN reaction_xref AS xref USING (rhea_id)
            ORDER BY
                group_id NULLS FIRST, input_id, rhea_id,
                reference_database, reference_id
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_CROSS_REFERENCE,
            output_schema=_public_schema(_SCHEMA_CROSS_REFERENCE, grouped=plan.grouped),
        )

    def uniprot_mappings(self) -> pl.LazyFrame:
        """Return one normalized row per selected reaction-to-UniProt edge.

        Examples:
            >>> lf = selection.uniprot_mappings()  # doctest: +SKIP
            >>> lf.select("rhea_id", "uniprot_id").collect()  # doctest: +SKIP
            shape: (..., 2)
        """

        plan = self._lazy_plan()
        _require_tables(
            plan,
            {"reaction_uniprot"},
            operation="read UniProt reaction mappings",
        )
        query = """
            WITH selected_reaction AS (
                SELECT DISTINCT rhea_id, master_id, direction
                FROM _selected_reaction
            )
            SELECT DISTINCT
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                mapping.uniprot_id AS uniprot_id,
                mapping.uniprot_section AS uniprot_section
            FROM selected_reaction AS selected
            JOIN reaction_uniprot AS mapping USING (rhea_id)
            ORDER BY rhea_id, uniprot_id, uniprot_section
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_UNIPROT_MAPPING,
            output_schema=_SCHEMA_UNIPROT_MAPPING,
        )

    def uniprot_neighborhoods(self) -> pl.LazyFrame:
        """Return one reaction row with nested input and UniProt lists.

        The two lists are assembled from independent ordered streams. The
        implementation never creates the input-by-UniProt Cartesian product.

        Examples:
            >>> lf = selection.uniprot_neighborhoods()  # doctest: +SKIP
            >>> by_input = lf.explode("inputs").unnest("inputs")  # doctest: +SKIP
            >>> by_input.collect()  # doctest: +SKIP
            shape: (..., ...)
        """

        plan = self._lazy_plan()
        _require_tables(
            plan,
            {"reaction_uniprot"},
            operation="read UniProt neighborhoods",
        )
        schema = _nested_neighborhood_schema(grouped=plan.grouped)
        return register_replayable_source(
            schema=schema,
            batches=lambda batch_size: _iter_neighborhood_batches(
                plan,
                batch_size=batch_size,
            ),
        )

    def publications(self) -> pl.LazyFrame:
        """Return PubMed references owned by selected reactions lazily.

        Examples:
            >>> lf = selection.publications()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        plan = self._lazy_plan()
        _require_tables(plan, {"reaction_publication"}, operation="read publications")
        query = """
            SELECT
                selected.group_id AS group_id,
                selected.input_id AS input_id,
                selected.input_namespace AS input_namespace,
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                publication.pubmed_id AS pubmed_id
            FROM _selected_reaction AS selected
            JOIN reaction_publication AS publication USING (rhea_id)
            ORDER BY group_id NULLS FIRST, input_id, rhea_id, pubmed_id
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_PUBLICATION,
            output_schema=_public_schema(_SCHEMA_PUBLICATION, grouped=plan.grouped),
        )

    def relationships(self) -> pl.LazyFrame:
        """Return hierarchy edges touching selected master reactions lazily.

        Examples:
            >>> lf = selection.relationships()  # doctest: +SKIP
            >>> lf.collect().head(1)  # doctest: +SKIP
            shape: (1, ...)
        """

        plan = self._lazy_plan()
        _require_tables(
            plan,
            {"reaction_relationship"},
            operation="read reaction relationships",
        )
        query = """
            SELECT DISTINCT
                selected.group_id AS group_id,
                selected.input_id AS input_id,
                selected.input_namespace AS input_namespace,
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                relation.from_reaction_id AS from_reaction_id,
                relation.to_reaction_id AS to_reaction_id,
                relation.relation_type AS relation_type
            FROM _selected_reaction AS selected
            JOIN reaction_relationship AS relation
              ON relation.from_reaction_id = selected.master_id
              OR relation.to_reaction_id = selected.master_id
            ORDER BY
                group_id NULLS FIRST, input_id, rhea_id,
                from_reaction_id, to_reaction_id, relation_type
        """
        return _relation_frame(
            plan,
            query=query,
            raw_schema=_SCHEMA_RELATIONSHIP,
            output_schema=_public_schema(_SCHEMA_RELATIONSHIP, grouped=plan.grouped),
        )

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return selected input IDs with no reaction match, lazily.

        Examples:
            >>> lf = selection.unmatched_ids()  # doctest: +SKIP
            >>> lf.collect()  # doctest: +SKIP
            shape: (..., ...)
        """

        if self._is_grouped:
            input_frame = pl.DataFrame(
                self._group_membership,
                schema={"group_id": pl.String, "input_id": pl.String},
                orient="row",
            ).lazy()
            matched = self.matches().select("group_id", "input_id")
            return (
                input_frame.join(matched, on=["group_id", "input_id"], how="anti")
                .with_columns(pl.lit(self.namespace).alias("input_namespace"))
                .select(list(_SCHEMA_UNMATCHED))
                .sort(list(_SCHEMA_UNMATCHED))
            )
        input_frame = pl.DataFrame(
            [(row.input_id, self.namespace) for row in self._input_rows],
            schema=_SCHEMA_UNIQUE_UNMATCHED,
            orient="row",
        ).lazy()
        matched = self.matches().select("input_id")
        return input_frame.join(matched, on="input_id", how="anti").sort(
            list(_SCHEMA_UNIQUE_UNMATCHED)
        )

    def _lazy_plan(self) -> _RheaLazyPlan:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]
        return _RheaLazyPlan(
            path=publication.path,
            identity=publication.identity,
            tables=publication.tables,
            views=publication.views,
            namespace=self.namespace,
            include_obsolete=self.include_obsolete,
            input_rows=self._input_rows,
            group_membership=self._group_membership,
            grouped=self._is_grouped,
        )

    def _eager_reactions(self) -> pl.DataFrame:
        """Extract selected reaction facts, including optional reaction SMILES.

        Examples:
            Read exact direction alongside the reaction equation:

            >>> selection.reactions().select(  # doctest: +SKIP
            ...     "rhea_id", "direction", "equation"
            ... ).collect_schema().names()
            ['rhea_id', 'direction', 'equation']
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
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                reaction.rhea_id AS rhea_id,
                reaction.master_id AS master_id,
                reaction.direction AS direction,
                reaction.accession AS accession,
                reaction.equation AS equation,
                reaction.equation_html AS equation_html,
                reaction.status AS status,
                reaction.is_balanced AS is_balanced,
                reaction.is_transport AS is_transport,
                reaction.comment AS comment,
                reaction.is_obsolete AS is_obsolete,
                {smiles_expression} AS reaction_smiles
            FROM _selected_reaction AS selected
            JOIN reaction USING (rhea_id)
            {smiles_join}
            ORDER BY group_id, input_id, rhea_id
        """
        return self._query_selected(query, schema=_SCHEMA_REACTION)

    def _eager_participants(self) -> pl.DataFrame:
        """Extract participants with compound facts and direction-aware roles.

        Examples:
            Read retained sides and nullable direction-specific roles:

            >>> selection.participants().select(  # doctest: +SKIP
            ...     "side", "directional_role", "chebi_id"
            ... ).collect_schema().names()
            ['side', 'directional_role', 'chebi_id']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_participant", "compound"},
            operation="extract reaction participants",
        )
        if "reaction_participant_direction" not in publication.views:
            raise CapabilityError(
                "Rhea publication lacks view required to extract reaction "
                "participants: reaction_participant_direction"
            )
        query = """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                participant.rhea_id AS rhea_id,
                participant.master_id AS master_id,
                participant.direction AS direction,
                participant.participant_id AS participant_id,
                participant.compound_id AS compound_id,
                participant.side AS side,
                participant.directional_role AS directional_role,
                participant.coefficient_text AS coefficient_text,
                participant.coefficient_numeric AS coefficient_numeric,
                participant.location AS location,
                compound.rhea_compound_id AS rhea_compound_id,
                compound.public_accession AS compound_accession,
                compound.compound_type AS compound_type,
                compound.name AS compound_name,
                compound.formula AS formula,
                compound.charge_text AS charge_text,
                compound.charge_numeric AS charge_numeric,
                compound.chebi_id AS chebi_id,
                compound.underlying_chebi_id AS underlying_chebi_id,
                compound.polymerization_index AS polymerization_index
            FROM _selected_reaction AS selected
            JOIN reaction_participant_direction AS participant USING (rhea_id)
            LEFT JOIN compound USING (compound_id)
            ORDER BY
                group_id, input_id, rhea_id, side, participant_id
        """
        return self._query_selected(query, schema=_SCHEMA_PARTICIPANT)

    def _eager_cross_references(self) -> pl.DataFrame:
        """Extract Rhea-owned external and UniProt reaction mappings.

        Examples:
            Materialize the normalized reference relation:

            >>> selection.cross_references().select(  # doctest: +SKIP
            ...     "reference_database", "reference_id"
            ... ).collect_schema().names()
            ['reference_database', 'reference_id']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        has_xref = "reaction_xref" in publication.tables
        has_uniprot = "reaction_uniprot" in publication.tables
        if not has_xref and not has_uniprot:
            raise CapabilityError(
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
                group_id AS "group_id",
                input_id AS "input_id",
                input_namespace AS "input_namespace",
                rhea_id AS rhea_id,
                master_id AS master_id,
                direction AS direction,
                reference_database AS reference_database,
                reference_id AS reference_id,
                uniprot_section AS uniprot_section
            FROM ({" UNION ALL ".join(branches)})
            ORDER BY
                group_id, input_id, rhea_id,
                reference_database, reference_id
        """
        return self._query_selected(query, schema=_SCHEMA_CROSS_REFERENCE)

    def _eager_publications(self) -> pl.DataFrame:
        """Extract PubMed references owned by selected Rhea reactions.

        Examples:
            Read PubMed IDs while retaining reaction lineage:

            >>> selection.publications().select(  # doctest: +SKIP
            ...     "rhea_id", "pubmed_id"
            ... ).collect_schema().names()
            ['rhea_id', 'pubmed_id']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_publication"},
            operation="extract reaction publications",
        )
        query = """
            SELECT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                publication.pubmed_id AS pubmed_id
            FROM _selected_reaction AS selected
            JOIN reaction_publication AS publication USING (rhea_id)
            ORDER BY group_id, input_id, rhea_id, pubmed_id
        """
        return self._query_selected(query, schema=_SCHEMA_PUBLICATION)

    def _eager_relationships(self) -> pl.DataFrame:
        """Extract hierarchy edges touching selected master reactions.

        Examples:
            Preserve each directed Rhea hierarchy edge:

            >>> selection.relationships().select(  # doctest: +SKIP
            ...     "from_reaction_id", "to_reaction_id", "relation_type"
            ... ).collect_schema().names()
            ['from_reaction_id', 'to_reaction_id', 'relation_type']
        """
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        _require_tables(
            publication,
            {"reaction_relationship"},
            operation="extract reaction relationships",
        )
        query = """
            SELECT DISTINCT
                selected.group_id AS "group_id",
                selected.input_id AS "input_id",
                selected.input_namespace AS "input_namespace",
                selected.rhea_id AS rhea_id,
                selected.master_id AS master_id,
                selected.direction AS direction,
                relation.from_reaction_id AS from_reaction_id,
                relation.to_reaction_id AS to_reaction_id,
                relation.relation_type AS relation_type
            FROM _selected_reaction AS selected
            JOIN reaction_relationship AS relation
              ON relation.from_reaction_id = selected.master_id
              OR relation.to_reaction_id = selected.master_id
            ORDER BY
                group_id, input_id, rhea_id,
                from_reaction_id, to_reaction_id, relation_type
        """
        return self._query_selected(query, schema=_SCHEMA_RELATIONSHIP)

    def _eager_unmatched_ids(self) -> pl.DataFrame:
        """Return normalized input identifiers with no matching Rhea reaction.

        Examples:
            Audit identifiers that did not resolve:

            >>> selection.unmatched_ids().collect_schema().names()  # doctest: +SKIP
            ['input_id', 'input_namespace']
        """
        matched_input_ids = set(self._eager_matches()["input_id"].to_list())
        rows = [
            (row.input_id, self.namespace)
            for row in self._input_rows
            if row.input_id not in matched_input_ids
        ]
        frame = pl.DataFrame(
            rows,
            schema=_SCHEMA_UNIQUE_UNMATCHED,
            orient="row",
        )
        frame = self._expand_groups(frame)
        expected = (
            list(_SCHEMA_UNMATCHED)
            if self._is_grouped
            else list(_SCHEMA_UNIQUE_UNMATCHED)
        )
        return frame.select(expected).sort(expected)

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
                input.input_id AS "input_id",
                '{self.namespace}' AS "input_namespace",
                reaction.rhea_id AS rhea_id,
                reaction.master_id AS master_id,
                reaction.direction AS direction
            FROM _input_id AS input
            {join}
            WHERE true {obsolete_filter}
            ORDER BY input_id, rhea_id
        """
        with _connect(publication) as connection:
            _create_input_table(connection, self._input_rows)
            frame = _fetch_frame(connection, query, schema=_SCHEMA_UNIQUE_MATCH)
        frame = self._expand_groups(frame)
        expected = (
            list(_SCHEMA_MATCH) if self._is_grouped else list(_SCHEMA_UNIQUE_MATCH)
        )
        return frame.select(expected).sort(expected)

    def _expand_groups(self, frame: pl.DataFrame) -> pl.DataFrame:
        if not self._is_grouped:
            return frame
        membership = pl.DataFrame(
            self._group_membership,
            schema={"group_id": pl.String, "input_id": pl.String},
            orient="row",
        )
        return membership.join(frame, on="input_id", how="inner").select(
            "group_id",
            *frame.columns,
        )

    def _query_selected(
        self,
        query: str,
        *,
        schema: SchemaDict,
    ) -> pl.DataFrame:
        publication = self.database._require_publication()  # pyright: ignore[reportPrivateUsage]  # sibling query boundary
        matches = self._eager_matches()
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
            expected.remove("group_id")
            frame = frame.drop("group_id")
        return frame.select(expected)


def open_rhea_publication(path: str | Path) -> _RheaPublication:
    """Validate and describe one bioextract Rhea DuckDB publication."""
    publication_path = Path(path)
    if not publication_path.exists():
        raise FileNotFoundError(publication_path)
    if not publication_path.is_file():
        raise ValueError(f"Rhea DuckDB path is not a file: {publication_path}")
    identity_before = _file_identity(publication_path)

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
        validate_duckdb_metadata_v1(connection, metadata)
        if metadata.get("bioextract.source_schema_profile") != SOURCE_SCHEMA_PROFILE:
            raise ValueError("Unsupported Rhea source schema profile")
        resource_schema_key = "bioextract.resource_schema_version"
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
    identity_after = _file_identity(publication_path)
    if identity_after != identity_before:
        raise IntegrityError("Rhea publication changed during validation")
    return _RheaPublication(
        path=publication_path,
        identity=identity_after,
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
        _InputRow(input_id=input_id, lookup_value=lookup_value)
        for input_id, lookup_value in _normalize_ids(ids, normalized_namespace)
    )
    return RheaReactionSelection(
        database=database,
        namespace=normalized_namespace,
        include_obsolete=bool(include_obsolete),
        _input_rows=rows,
        _group_ids=(),
        _group_membership=(),
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
    normalized_group_ids = [str(group_id).strip() for group_id in ids_by_group]
    validate_group_ids(normalized_group_ids)
    unique_rows: set[_InputRow] = set()
    membership: set[tuple[str, str]] = set()
    for raw_group_id, ids in ids_by_group.items():
        group_id = str(raw_group_id).strip()
        for input_id, lookup_value in _normalize_ids(ids, normalized_namespace):
            unique_rows.add(_InputRow(input_id=input_id, lookup_value=lookup_value))
            membership.add((group_id, input_id))
    rows = sorted(
        unique_rows,
        key=lambda row: (row.input_id, row.lookup_value),
    )
    return RheaReactionSelection(
        database=database,
        namespace=normalized_namespace,
        include_obsolete=bool(include_obsolete),
        _input_rows=tuple(rows),
        _group_ids=tuple(sorted(normalized_group_ids)),
        _group_membership=tuple(sorted(membership)),
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
    publication: _RheaPublication | _RheaLazyPlan,
    required: set[str],
    *,
    operation: str,
) -> None:
    missing = sorted(required - set(publication.tables))
    if missing:
        raise CapabilityError(
            f"Rhea publication cannot {operation}; missing relations: {missing}"
        )


def _public_schema(schema: SchemaDict, *, grouped: bool) -> SchemaDict:
    if grouped:
        return dict(schema)
    return {name: dtype for name, dtype in schema.items() if name != "group_id"}


def _relation_frame(
    plan: _RheaLazyPlan,
    *,
    query: str,
    raw_schema: SchemaDict,
    output_schema: SchemaDict,
    prepare_selected: bool = True,
) -> pl.LazyFrame:
    return register_replayable_source(
        schema=output_schema,
        batches=lambda batch_size: _iter_query_batches(
            plan,
            query=query,
            raw_schema=raw_schema,
            batch_size=batch_size,
            prepare_selected=prepare_selected,
        ),
    )


def _matches_query(plan: _RheaLazyPlan) -> str:
    obsolete_filter = "" if plan.include_obsolete else "AND NOT reaction.is_obsolete"
    if plan.namespace == "rhea":
        join = """
            JOIN reaction
              ON CAST(reaction.rhea_id AS VARCHAR) = input.lookup_value
        """
    elif plan.namespace == "chebi":
        join = """
            JOIN compound
              ON compound.chebi_id = input.lookup_value
            JOIN reaction_participant AS participant
              ON participant.compound_id = compound.compound_id
            JOIN reaction
              ON reaction.master_id = participant.master_id
        """
    elif plan.namespace == "uniprot":
        join = """
            JOIN reaction_uniprot AS mapping
              ON mapping.uniprot_id = input.lookup_value
            JOIN reaction USING (rhea_id)
        """
    else:
        database_name = _EXTERNAL_DATABASE_BY_NAMESPACE[plan.namespace]
        join = f"""
            JOIN reaction_xref AS mapping
              ON mapping.external_id = input.lookup_value
             AND mapping.database_name = '{database_name}'
            JOIN reaction USING (rhea_id)
        """
    group_select = (
        "membership.group_id AS group_id"
        if plan.grouped
        else "CAST(NULL AS VARCHAR) AS group_id"
    )
    group_join = (
        "LEFT JOIN _group_membership AS membership USING (input_id)"
        if plan.grouped
        else ""
    )
    return f"""
        SELECT DISTINCT
            {group_select},
            input.input_id AS input_id,
            '{plan.namespace}' AS input_namespace,
            reaction.rhea_id AS rhea_id,
            reaction.master_id AS master_id,
            reaction.direction AS direction
        FROM _input_id AS input
        {group_join}
        {join}
        WHERE true {obsolete_filter}
        ORDER BY group_id NULLS FIRST, input_id, rhea_id
    """


def _iter_query_batches(
    plan: _RheaLazyPlan,
    *,
    query: str,
    raw_schema: SchemaDict,
    batch_size: int,
    prepare_selected: bool,
) -> Iterator[pl.DataFrame]:
    connection = _connect_checked(plan)
    try:
        _create_input_table(connection, plan.input_rows)
        _create_group_membership_table(connection, plan)
        if prepare_selected:
            _prepare_selected_reaction(connection, plan)
        result = connection.execute(query)
        reader = _arrow_reader(result, batch_size)
        try:
            for record_batch in reader:
                frame: pl.DataFrame = pl.from_arrow(record_batch)  # type: ignore[reportUnknownMemberType]
                frame = frame.cast(raw_schema, strict=False)  # type: ignore[reportArgumentType]
                if not plan.grouped and "group_id" in frame.columns:
                    frame = frame.drop("group_id")
                yield frame
        finally:
            close = getattr(reader, "close", None)
            if close is not None:
                close()
    finally:
        connection.close()


def _iter_neighborhood_batches(
    plan: _RheaLazyPlan,
    *,
    batch_size: int,
) -> Iterator[pl.DataFrame]:
    connection = _connect_checked(plan)
    mapping_connection = _connect_checked(plan)
    output_schema = _nested_neighborhood_schema(grouped=plan.grouped)
    try:
        _create_input_table(connection, plan.input_rows)
        _create_group_membership_table(connection, plan)
        _prepare_selected_reaction(connection, plan)
        _create_input_table(mapping_connection, plan.input_rows)
        _create_group_membership_table(mapping_connection, plan)
        _prepare_selected_reaction(mapping_connection, plan)
        matches_reader = _arrow_reader(
            connection.execute(
                """
                SELECT DISTINCT
                    rhea_id, master_id, direction, group_id,
                    input_id, input_namespace
                FROM _selected_reaction
                ORDER BY rhea_id, group_id NULLS FIRST, input_id
                """
            ),
            batch_size,
        )
        uniprot_reader = _arrow_reader(
            mapping_connection.execute(
                """
                SELECT DISTINCT
                    rhea_id, uniprot_id, uniprot_section
                FROM reaction_uniprot
                JOIN (
                    SELECT DISTINCT rhea_id FROM _selected_reaction
                ) AS selected USING (rhea_id)
                ORDER BY rhea_id, uniprot_id, uniprot_section
                """
            ),
            batch_size,
        )
        try:
            match_rows = _iter_arrow_rows(matches_reader)
            uniprot_rows = _iter_arrow_rows(uniprot_reader)
            next_match = next(match_rows, None)
            next_uniprot = next(uniprot_rows, None)
            output_rows: list[dict[str, object]] = []
            # The generic adapter may request a large Arrow batch. Nested
            # neighborhoods have deliberately high per-row cardinality, so
            # keep only a small number of completed reaction rows in memory.
            output_batch_size = min(batch_size, _NESTED_OUTPUT_BATCH_SIZE)
            while next_match is not None:
                rhea_id = int(str(next_match["rhea_id"]))
                master_id = next_match["master_id"]
                direction = next_match["direction"]
                inputs: list[dict[str, object]] = []
                while (
                    next_match is not None
                    and int(str(next_match["rhea_id"])) == rhea_id
                ):
                    input_row: dict[str, object] = {
                        "input_id": next_match["input_id"],
                        "input_namespace": next_match["input_namespace"],
                    }
                    if plan.grouped:
                        input_row = {
                            "group_id": next_match["group_id"],
                            **input_row,
                        }
                    inputs.append(input_row)
                    next_match = next(match_rows, None)

                uniprot_entries: list[dict[str, object]] = []
                while (
                    next_uniprot is not None
                    and int(str(next_uniprot["rhea_id"])) < rhea_id
                ):
                    next_uniprot = next(uniprot_rows, None)
                while (
                    next_uniprot is not None
                    and int(str(next_uniprot["rhea_id"])) == rhea_id
                ):
                    uniprot_entries.append(
                        {
                            "uniprot_id": next_uniprot["uniprot_id"],
                            "uniprot_section": next_uniprot["uniprot_section"],
                        }
                    )
                    next_uniprot = next(uniprot_rows, None)
                output_rows.append(
                    {
                        "rhea_id": rhea_id,
                        "master_id": master_id,
                        "direction": direction,
                        "inputs": inputs,
                        "uniprot_entries": uniprot_entries,
                    }
                )
                if len(output_rows) >= output_batch_size:
                    yield pl.DataFrame(output_rows, schema=output_schema, strict=False)
                    output_rows = []
            if output_rows:
                yield pl.DataFrame(output_rows, schema=output_schema, strict=False)
        finally:
            for reader in (matches_reader, uniprot_reader):
                close = getattr(reader, "close", None)
                if close is not None:
                    close()
    finally:
        mapping_connection.close()
        connection.close()


def _connect_checked(plan: _RheaLazyPlan) -> duckdb.DuckDBPyConnection:
    if _file_identity(plan.path) != plan.identity:
        raise IntegrityError(
            "Rhea publication was replaced; reopen it with from_duckdb()"
        )
    try:
        connection = duckdb.connect(str(plan.path), read_only=True)
    except duckdb.Error as error:
        raise IntegrityError(
            "Rhea publication became unavailable; reopen it with from_duckdb()"
        ) from error
    if _file_identity(plan.path) != plan.identity:
        connection.close()
        raise IntegrityError(
            "Rhea publication was replaced; reopen it with from_duckdb()"
        )
    return connection


def _arrow_reader(result: Any, batch_size: int) -> Any:
    """Use the current DuckDB Arrow reader and retain an older fallback."""

    to_arrow_reader = getattr(result, "to_arrow_reader", None)
    if to_arrow_reader is not None:
        return to_arrow_reader(batch_size)
    return result.fetch_record_batch(rows_per_batch=batch_size)


def _iter_arrow_rows(reader: Any) -> Iterator[dict[str, object]]:
    for record_batch in reader:
        yield from record_batch.to_pylist()


def _create_group_membership_table(
    connection: duckdb.DuckDBPyConnection,
    plan: _RheaLazyPlan,
) -> None:
    if not plan.grouped:
        return
    connection.execute(
        """
        CREATE TEMP TABLE _group_membership (
            group_id VARCHAR NOT NULL,
            input_id VARCHAR NOT NULL
        )
        """
    )
    if plan.group_membership:
        connection.executemany(
            "INSERT INTO _group_membership VALUES (?, ?)",
            list(plan.group_membership),
        )


def _prepare_selected_reaction(
    connection: duckdb.DuckDBPyConnection,
    plan: _RheaLazyPlan,
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
    reader = _arrow_reader(connection.execute(_matches_query(plan)), 100_000)
    group_by_input: dict[str, tuple[str, ...]] = {}
    if plan.grouped:
        groups: dict[str, list[str]] = {}
        for group_id, input_id in plan.group_membership:
            groups.setdefault(input_id, []).append(group_id)
        group_by_input = {
            input_id: tuple(group_ids) for input_id, group_ids in groups.items()
        }
    rows: list[tuple[object, ...]] = []
    try:
        for row in _iter_arrow_rows(reader):
            group_ids = (
                group_by_input.get(str(row["input_id"]), ())
                if plan.grouped
                else (None,)
            )
            for group_id in group_ids:
                rows.append(
                    (
                        group_id,
                        row["input_id"],
                        row["input_namespace"],
                        row["rhea_id"],
                        row["master_id"],
                        row["direction"],
                    )
                )
                if len(rows) >= 10_000:
                    connection.executemany(
                        "INSERT INTO _selected_reaction VALUES (?, ?, ?, ?, ?, ?)",
                        rows,
                    )
                    rows = []
        if rows:
            connection.executemany(
                "INSERT INTO _selected_reaction VALUES (?, ?, ?, ?, ?, ?)",
                rows,
            )
    finally:
        close = getattr(reader, "close", None)
        if close is not None:
            close()


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
            row["group_id"] if grouped else None,
            row["input_id"],
            row["input_namespace"],
            row["rhea_id"],
            row["master_id"],
            row["direction"],
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


def _file_identity(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
        int(stat.st_ctime_ns),
    )
