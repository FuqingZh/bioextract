from __future__ import annotations

import copy
import re
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from threading import Lock
from typing import Any

import duckdb
import polars as pl
from polars._typing import SchemaDict

from bioextract._lazy import (
    _RelationScanRequest,  # pyright: ignore[reportPrivateUsage]  # typed source boundary
    register_replayable_source,
)
from bioextract._shared import create_group_input_frames
from bioextract.errors import CapabilityError, IntegrityError

from .constant import (
    NAMESPACE_VALUES,
    SCHEMA_GENE_PATHWAY_VIA_KO,
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_MATCH,
    SCHEMA_ORGANISM,
    SCHEMA_UNMAPPED,
    TABLE_SCHEMAS,
    KEGGNamespace,
)
from .parse import aggregate_kos, aggregate_organism, parse_organism_list
from .source import MappingSnapshot, resolve_organism_work, source_capabilities

_UNIPROT_PIPE = re.compile(r"^[^|]+\|([^|]+)\|")


@dataclass(slots=True)
class _KeggMatchState:
    lock: Lock = field(default_factory=Lock)
    frame: pl.DataFrame | None = None


def relation(snapshot: MappingSnapshot, name: str) -> pl.LazyFrame:
    schema = _relation_schema(name)
    frozen = copy.copy(snapshot)
    if frozen.mode == "publication":
        return register_replayable_source(
            schema=schema,
            batches=lambda request: _iter_publication(
                frozen,
                name=name,
                schema=schema,
                request=request,
            ),
        )
    return register_replayable_source(
        schema=schema,
        batches=lambda request: _iter_source(frozen, name=name, request=request),
    )


def gene_pathways(snapshot: MappingSnapshot) -> pl.LazyFrame:
    _require_capabilities(snapshot, "gene_pathways", "gene_pathway")
    return relation(snapshot, "gene_annotation").select(
        "organism_code", "kegg_gene_id", "pathway_mappings"
    )


def ko_pathways(snapshot: MappingSnapshot) -> pl.LazyFrame:
    _require_capabilities(snapshot, "ko_pathways", "ko_pathway")
    return relation(snapshot, "ko_annotation")


def gene_pathways_via_ko(snapshot: MappingSnapshot) -> pl.LazyFrame:
    _require_capabilities(snapshot, "gene_pathways_via_ko", "gene_ko", "ko_pathway")
    genes = (
        relation(snapshot, "gene_annotation")
        .select("organism_code", "kegg_gene_id", "ko_mappings")
        .explode("ko_mappings")
        .drop_nulls("ko_mappings")
        .unnest("ko_mappings")
    )
    pathways = (
        relation(snapshot, "ko_annotation")
        .explode("pathway_mappings")
        .drop_nulls("pathway_mappings")
        .unnest("pathway_mappings")
    )
    return (
        genes.join(pathways, on="ko_id", how="inner")
        .group_by("organism_code", "kegg_gene_id")
        .agg(
            pl.struct(
                "ko_id",
                "kegg_pathway_id",
                "pathway_namespace",
                "pathway_map_id",
            )
            .unique()
            .sort()
            .alias("pathway_mappings")
        )
        .select(*SCHEMA_GENE_PATHWAY_VIA_KO)
    )


@dataclass(frozen=True, slots=True)
class KeggSelection:
    """Replayable caller-to-KEGG-gene selection over one mapping handle.

    Examples:
        >>> selection = db.select_ids(["P12345"], namespace="uniprot")  # doctest: +SKIP
        >>> selection.matches().collect()  # doctest: +SKIP
        shape: (..., 4)
    """

    snapshot: MappingSnapshot
    input_ids: pl.DataFrame
    namespace: KEGGNamespace
    group_membership: pl.DataFrame | None = None
    _match_state: _KeggMatchState = field(
        default_factory=_KeggMatchState,
        init=False,
        repr=False,
        compare=False,
    )

    @classmethod
    def _from_ids(
        cls,
        snapshot: MappingSnapshot,
        ids: Iterable[str],
        *,
        namespace: KEGGNamespace,
    ) -> KeggSelection:
        return cls(
            snapshot=copy.copy(snapshot),
            input_ids=_input_frame(ids, namespace=namespace),
            namespace=namespace,
        )

    @classmethod
    def _from_groups(
        cls,
        snapshot: MappingSnapshot,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: KEGGNamespace,
    ) -> KeggSelection:
        normalized = {
            group: [
                normalized_id
                for value in values
                if (normalized_id := _normalize_input(value, namespace=namespace))
            ]
            for group, values in ids_by_group.items()
        }
        frames = create_group_input_frames(
            normalized,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return cls(
            snapshot=copy.copy(snapshot),
            input_ids=frames.df_input_ids,
            namespace=namespace,
            group_membership=frames.df_group_membership,
        )

    @property
    def grouped(self) -> bool:
        """Report whether selection outputs retain caller group identity.

        Examples:
            >>> selection.grouped  # doctest: +SKIP
            False
        """
        return self.group_membership is not None

    def matches(self) -> pl.LazyFrame:
        """Return normalized input-to-gene lineage without annotation expansion.

        Examples:
            >>> selection.matches().collect()  # doctest: +SKIP
            shape: (..., 4)
        """
        capability = {
            "uniprot": "uniprot_conversion",
            "ncbi_gene": "ncbi_gene_conversion",
        }.get(self.namespace)
        if capability is not None:
            _require_capabilities(self.snapshot, "match selected IDs", capability)
        if self.snapshot.mode == "publication":
            snapshot = copy.copy(self)
            columns = ["input_id", "input_namespace", "organism_code", "kegg_gene_id"]
            if self.grouped:
                columns.insert(0, "group_id")
            schema: SchemaDict = {
                column: (
                    SCHEMA_GROUP_INPUT_IDS["group_id"]
                    if column == "group_id"
                    else SCHEMA_MATCH[column]
                )
                for column in columns
            }
            return register_replayable_source(
                schema=schema,
                batches=lambda request: _iter_publication_matches(
                    snapshot,
                    schema=schema,
                    request=request,
                ),
                is_pure=True,
            )
        genes = relation(self.snapshot, "gene_annotation")
        if self.namespace == "kegg_gene":
            lookup = genes.select("organism_code", "kegg_gene_id").unique()
            matched = (
                self.input_ids.lazy()
                .join(
                    lookup,
                    left_on="input_id",
                    right_on="kegg_gene_id",
                    how="inner",
                )
                .with_columns(pl.col("input_id").alias("kegg_gene_id"))
            )
        else:
            list_column, value_column = {
                "uniprot": ("uniprot_mappings", "uniprot_id"),
                "ncbi_gene": ("ncbi_gene_mappings", "ncbi_gene_id"),
            }[self.namespace]
            lookup = (
                genes.select("organism_code", "kegg_gene_id", list_column)
                .explode(list_column)
                .drop_nulls(list_column)
                .unnest(list_column)
            )
            matched = self.input_ids.lazy().join(
                lookup,
                left_on="input_id",
                right_on=value_column,
                how="inner",
            )
        matched = (
            matched.with_columns(pl.lit(self.namespace).alias("input_namespace"))
            .select(*SCHEMA_MATCH)
            .unique()
        )
        if self.group_membership is None:
            return matched
        return (
            self.group_membership.lazy()
            .join(matched, on="input_id", how="inner")
            .select("group_id", *SCHEMA_MATCH)
            .unique()
        )

    def gene_annotations(self) -> pl.LazyFrame:
        """Return one aggregate annotation row per matched composite gene.

        Examples:
            >>> selection.gene_annotations().collect()  # doctest: +SKIP
            shape: (..., 12)
        """
        return self._attach_inputs(relation(self.snapshot, "gene_annotation"))

    def gene_pathways(self) -> pl.LazyFrame:
        """Return direct pathway observations for each matched composite gene.

        Examples:
            >>> selection.gene_pathways().collect()  # doctest: +SKIP
            shape: (..., 4)
        """
        return self._attach_inputs(gene_pathways(self.snapshot))

    def gene_pathways_via_ko(self) -> pl.LazyFrame:
        """Return KO-mediated pathway evidence for each matched composite gene.

        Examples:
            >>> selection.gene_pathways_via_ko().collect()  # doctest: +SKIP
            shape: (..., 4)
        """
        return self._attach_inputs(gene_pathways_via_ko(self.snapshot))

    def unmatched_ids(self) -> pl.LazyFrame:
        """Return normalized IDs unmatched after the handle's organism scope.

        Examples:
            >>> selection.unmatched_ids().collect()  # doctest: +SKIP
            shape: (..., 1)
        """
        matched = self.matches().select("input_id").unique()
        unmatched = self.input_ids.lazy().join(matched, on="input_id", how="anti")
        if self.group_membership is None:
            return unmatched.select("input_id")
        return (
            self.group_membership.lazy()
            .join(unmatched, on="input_id", how="inner")
            .select("group_id", "input_id")
        )

    def _attach_inputs(self, annotation: pl.LazyFrame) -> pl.LazyFrame:
        fields = ["input_id", "input_namespace"]
        if self.grouped:
            fields.insert(0, "group_id")
        inputs = (
            self.matches()
            .group_by("organism_code", "kegg_gene_id")
            .agg(pl.struct(*fields).unique().sort().alias("inputs"))
        )
        columns = list(annotation.collect_schema().names())
        return annotation.join(
            inputs, on=["organism_code", "kegg_gene_id"], how="inner"
        ).select(*columns[:2], "inputs", *columns[2:])


def _iter_publication_matches(
    selection: KeggSelection,
    *,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    with selection._match_state.lock:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        if selection._match_state.frame is None:  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
            selection._match_state.frame = _query_publication_matches(selection)  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
        frame = selection._match_state.frame.clone()  # pyright: ignore[reportPrivateUsage]  # paired selection boundary
    requested = _select_requested(frame, request.columns)
    yield from _slices(requested, request.effective_batch_size)


def _query_publication_matches(selection: KeggSelection) -> pl.DataFrame:
    path = selection.snapshot.publication_path
    if path is None:
        raise IntegrityError("KEGG mapping publication path is missing")
    input_ids = selection.input_ids.get_column("input_id").to_list()
    schema = SCHEMA_MATCH
    if not input_ids:
        frame = pl.DataFrame(schema=schema)
    else:
        connection = duckdb.connect(str(path), read_only=True)
        try:
            connection.execute(
                "CREATE TEMP TABLE _selected_input(input_id VARCHAR PRIMARY KEY)"
            )
            connection.executemany(
                "INSERT INTO _selected_input VALUES (?)",
                [(str(value),) for value in input_ids],
            )
            parameters: list[str] = []
            scope = ""
            if selection.snapshot.organism_scope is not None:
                scope = (
                    " WHERE gene.organism_code IN ("
                    + ",".join("?" for _ in selection.snapshot.organism_scope)
                    + ")"
                )
                parameters.extend(selection.snapshot.organism_scope)
            if selection.namespace == "kegg_gene":
                query = (
                    """
                    SELECT DISTINCT
                        selected.input_id AS input_id,
                        'kegg_gene' AS input_namespace,
                        gene.organism_code,
                        gene.kegg_gene_id
                    FROM _selected_input AS selected
                    JOIN gene_annotation AS gene
                      ON gene.kegg_gene_id = selected.input_id
                """
                    + scope
                )
            else:
                list_column, value_column = {
                    "uniprot": ("uniprot_mappings", "uniprot_id"),
                    "ncbi_gene": ("ncbi_gene_mappings", "ncbi_gene_id"),
                }[selection.namespace]
                query = (
                    f"""
                    SELECT DISTINCT
                        selected.input_id AS input_id,
                        '{selection.namespace}' AS input_namespace,
                        gene.organism_code,
                        gene.kegg_gene_id
                    FROM gene_annotation AS gene
                    CROSS JOIN UNNEST(gene."{list_column}") AS mapping(item)
                    JOIN _selected_input AS selected
                      ON mapping.item."{value_column}" = selected.input_id
                """
                    + scope
                )
            rows = connection.execute(query, parameters).fetchall()
            frame = pl.DataFrame(rows, schema=schema, orient="row")
        finally:
            connection.close()
    if selection.group_membership is None:
        return frame.sort(list(schema))
    return (
        selection.group_membership.lazy()
        .join(frame.lazy(), on="input_id", how="inner")
        .select("group_id", *schema)
        .unique()
        .collect()
        .sort(["group_id", *schema])
    )


def _iter_source(
    snapshot: MappingSnapshot,
    *,
    name: str,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    batch_size = request.effective_batch_size
    metadata, _, _ = parse_organism_list(snapshot.organism_list)
    work = resolve_organism_work(snapshot, validate_role_files=name == "organism")
    if name == "organism":
        rows: list[dict[str, object]] = [
            {
                "organism_code": code,
                "genome_id": metadata.get(code, {}).get("genome_id"),
                "organism_name": metadata.get(code, {}).get("organism_name"),
                "taxonomy_lineage": (
                    metadata.get(code, {}).get("taxonomy_lineage")
                    if code in metadata
                    else ([] if snapshot.organism_list is not None else None)
                ),
            }
            for code, _ in work
        ]
        frame = pl.DataFrame(rows, schema=SCHEMA_ORGANISM)
        yield from _slices(_select_requested(frame, request.columns), batch_size)
        return
    if name == "gene_annotation":
        for code, roles in work:
            aggregate = aggregate_organism(
                code,
                roles,
                organism_metadata=metadata.get(code),
                organism_list_available=snapshot.organism_list is not None,
            )
            yield from _slices(
                _select_requested(aggregate.genes, request.columns), batch_size
            )
        return
    if name == "ko_annotation":
        ko_ids: set[str] = set()
        for code, roles in work:
            aggregate = aggregate_organism(
                code,
                roles,
                organism_metadata=metadata.get(code),
                organism_list_available=snapshot.organism_list is not None,
            )
            ko_ids.update(aggregate.ko_ids)
        frame, _, _ = aggregate_kos(ko_ids, snapshot.ko_pathway)
        yield from _slices(_select_requested(frame, request.columns), batch_size)
        return
    raise AssertionError(name)


def _iter_publication(
    snapshot: MappingSnapshot,
    *,
    name: str,
    schema: SchemaDict,
    request: _RelationScanRequest,
) -> Iterator[pl.DataFrame]:
    batch_size = request.effective_batch_size
    path = snapshot.publication_path
    if path is None:
        raise IntegrityError("KEGG mapping publication path is missing")
    try:
        connection = duckdb.connect(str(path), read_only=True)
    except duckdb.Error as error:
        raise IntegrityError(
            f"Cannot reopen KEGG mapping publication: {path}"
        ) from error
    try:
        parameters: list[str] = []
        where = ""
        if snapshot.organism_scope is not None and name in {
            "organism",
            "gene_annotation",
        }:
            placeholders = ", ".join("?" for _ in snapshot.organism_scope)
            where = f" WHERE organism_code IN ({placeholders})"
            parameters.extend(snapshot.organism_scope)
        columns = [
            column for column in (request.columns or tuple(schema)) if column in schema
        ]
        if not columns:
            columns = list(schema)
        projection = ", ".join(f'"{column}"' for column in columns)
        result = connection.execute(
            f'SELECT {projection} FROM "{name}"{where}',
            parameters,
        )
        reader = _arrow_reader(result, batch_size)
        for batch in reader:
            frame = pl.from_arrow(batch)  # pyright: ignore[reportUnknownMemberType]
            if not isinstance(frame, pl.DataFrame):
                raise IntegrityError("KEGG mapping query returned a non-tabular batch")
            yield frame.cast(  # type: ignore[reportArgumentType]
                {column: schema[column] for column in frame.columns},
                strict=False,
            )
    finally:
        connection.close()


def _relation_schema(name: str) -> SchemaDict:
    try:
        return dict(TABLE_SCHEMAS[name])
    except KeyError as error:
        raise ValueError(f"Unknown KEGG mapping relation: {name}") from error


def _arrow_reader(result: Any, batch_size: int) -> Any:
    to_arrow_reader = getattr(result, "to_arrow_reader", None)
    if to_arrow_reader is not None:
        return to_arrow_reader(batch_size)
    return result.fetch_record_batch(rows_per_batch=batch_size)  # pyright: ignore[reportUnknownMemberType]


def _slices(frame: pl.DataFrame, batch_size: int) -> Iterator[pl.DataFrame]:
    for offset in range(0, frame.height, batch_size):
        yield frame.slice(offset, batch_size)


def _select_requested(
    frame: pl.DataFrame,
    columns: tuple[str, ...] | None,
) -> pl.DataFrame:
    if columns is None:
        return frame
    selected = [column for column in columns if column in frame.columns]
    return frame.select(selected) if selected else frame


def _require_capabilities(
    snapshot: MappingSnapshot, operation: str, *required: str
) -> None:
    capabilities = source_capabilities(snapshot)
    missing = [name for name in required if not capabilities.get(name, False)]
    if missing:
        raise CapabilityError(
            f"KEGG mapping cannot {operation}; unavailable capabilities: {missing}"
        )


def validate_namespace(namespace: str) -> None:
    if namespace not in NAMESPACE_VALUES:
        raise ValueError(
            "namespace must be one of: "
            f"{', '.join(NAMESPACE_VALUES)}; got {namespace!r}"
        )


def _input_frame(ids: Iterable[str], *, namespace: KEGGNamespace) -> pl.DataFrame:
    values = [
        normalized
        for value in ids
        if (normalized := _normalize_input(value, namespace=namespace))
    ]
    if not values:
        return pl.DataFrame(schema=SCHEMA_UNMAPPED)
    return (
        pl.DataFrame({"input_id": values}, schema=SCHEMA_UNMAPPED)
        .unique()
        .sort("input_id")
    )


def _normalize_input(value: Any, *, namespace: KEGGNamespace) -> str:
    validate_namespace(namespace)
    normalized = str(value).strip()
    if not normalized:
        return ""
    prefixes = ("up:", "ncbi-geneid:", "ko:", "path:")
    if namespace == "uniprot":
        if match := _UNIPROT_PIPE.match(normalized):
            normalized = match.group(1).strip()
        elif normalized.startswith("up:"):
            normalized = normalized[3:]
    elif namespace == "ncbi_gene" and normalized.startswith("ncbi-geneid:"):
        normalized = normalized[len("ncbi-geneid:") :]
    elif namespace == "kegg_gene":
        prefix, separator, suffix = normalized.partition(":")
        if not separator or not suffix or re.fullmatch(r"[a-z]{3,4}", prefix) is None:
            raise ValueError(
                f"KEGG gene input requires a full organism prefix: {value!r}"
            )
    if any(normalized.startswith(prefix) for prefix in prefixes):
        raise ValueError(
            f"Input identifier prefix does not match namespace {namespace!r}: {value!r}"
        )
    return normalized
