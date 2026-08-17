"""Private bridges from replayable batch readers to Polars lazy sources.

The public API exposes the native :class:`polars.LazyFrame`; this module owns
only the execution boundary needed to make a resource-backed reader replayable
and to close its resources when Polars stops consuming it.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass

import polars as pl
from polars._typing import SchemaDict
from polars.exceptions import ComputeError

_DEFAULT_BATCH_SIZE = 100_000


@dataclass(frozen=True, slots=True)
class _RelationScanRequest:
    """Execution hints negotiated between Polars and a resource producer."""

    columns: tuple[str, ...] | None
    predicate: pl.Expr | None
    n_rows: int | None
    batch_size: int | None

    @property
    def effective_batch_size(self) -> int:
        """Return a positive batch hint for producers."""
        return self.batch_size or _DEFAULT_BATCH_SIZE


BatchFactory = Callable[[_RelationScanRequest], Iterator[pl.DataFrame]]
FrameFactory = Callable[[], pl.DataFrame]


def register_replayable_source(
    *,
    schema: SchemaDict,
    batches: BatchFactory,
    validate_schema: bool = True,
    is_pure: bool = False,
) -> pl.LazyFrame:
    """Create a native lazy frame backed by a replayable batch factory.

    ``batches`` is called once per Polars execution with a
    :class:`_RelationScanRequest` and must open all source resources inside the
    returned iterator. The iterator is closed when execution succeeds, fails,
    or is stopped early. Projection, predicates, and explicit row limits are
    applied to each emitted batch so the source remains correct when Polars
    pushes those operations into the IO plugin. Resource producers may use the
    request to push safe constraints into their native scanner or SQL query.

    This function intentionally stays private. Polars' IO-plugin API is marked
    unstable; resource modules depend on this small boundary rather than
    exposing plugin details to callers.
    """

    frozen_schema = dict(schema)

    def io_source(
        with_columns: list[str] | None,
        predicate: pl.Expr | None,
        n_rows: int | None,
        batch_size: int | None,
    ) -> Iterator[pl.DataFrame]:
        request = _RelationScanRequest(
            columns=_required_columns(with_columns, predicate),
            predicate=predicate,
            n_rows=n_rows,
            batch_size=batch_size,
        )
        if n_rows == 0:
            return
        source_iterator = batches(request)
        yielded_rows = 0
        try:
            for frame in source_iterator:
                if predicate is not None:
                    frame = frame.filter(predicate)
                if with_columns is not None:
                    frame = frame.select(with_columns)
                if n_rows is not None:
                    remaining = n_rows - yielded_rows
                    if remaining <= 0:
                        return
                    frame = frame.head(remaining)
                if frame.is_empty():
                    continue
                yielded_rows += frame.height
                yield frame
                if n_rows is not None and yielded_rows >= n_rows:
                    return
        finally:
            close = getattr(source_iterator, "close", None)
            if close is not None:
                close()

    # register_io_source is intentionally isolated here because Polars marks
    # this API unstable. The runtime dependency is checked at import time by
    # the repository's compatibility tests.
    return pl.io.plugins.register_io_source(  # type: ignore[reportAttributeAccessIssue, reportUnknownVariableType]
        io_source,
        schema=frozen_schema,
        validate_schema=validate_schema,
        is_pure=is_pure,
    )


def register_materialized_source(
    *,
    columns: Iterable[str] | None = None,
    schema: SchemaDict | None = None,
    frame: FrameFactory,
) -> pl.LazyFrame:
    """Defer a parser-backed frame factory until Polars executes the plan.

    This explicit fallback is for small or inherently global relations whose
    algorithm already owns a complete domain transformation but does not
    expose a batch reader. Column names are frozen for plan construction;
    concrete dtypes are supplied by the emitted frame. High-cardinality
    relations must use :func:`register_replayable_source` with an explicit
    typed schema instead.
    """

    if schema is None:
        if columns is None:
            raise ValueError("columns or schema must be provided")
        frozen_schema = dict.fromkeys(columns, pl.Unknown)
    else:
        if columns is not None:
            raise ValueError("columns and schema are mutually exclusive")
        frozen_schema = dict(schema)

    def batches(request: _RelationScanRequest) -> Iterator[pl.DataFrame]:
        del request
        yield frame()

    return register_replayable_source(
        schema=frozen_schema,
        batches=batches,
        validate_schema=False,
    )


def _required_columns(
    with_columns: list[str] | None,
    predicate: pl.Expr | None,
) -> tuple[str, ...] | None:
    """Return producer columns needed for final projection and filtering."""
    if with_columns is None:
        return None
    columns = list(dict.fromkeys(with_columns))
    if predicate is None:
        return tuple(columns)
    try:
        predicate_columns = predicate.meta.root_names()
    except (AttributeError, ComputeError, TypeError):
        # A producer that cannot reason about the expression must emit its
        # complete schema; the adapter still applies the predicate correctly.
        return None
    for name in predicate_columns:
        if name not in columns:
            columns.append(name)
    return tuple(columns)
