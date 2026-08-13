"""Private bridges from replayable batch readers to Polars lazy sources.

The public API exposes the native :class:`polars.LazyFrame`; this module owns
only the execution boundary needed to make a resource-backed reader replayable
and to close its resources when Polars stops consuming it.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator

import polars as pl
from polars._typing import SchemaDict

_DEFAULT_BATCH_SIZE = 100_000
BatchFactory = Callable[[int], Iterator[pl.DataFrame]]
FrameFactory = Callable[[], pl.DataFrame]


def register_replayable_source(
    *,
    schema: SchemaDict,
    batches: BatchFactory,
    validate_schema: bool = True,
) -> pl.LazyFrame:
    """Create a native lazy frame backed by a replayable batch factory.

    ``batches`` is called once per Polars execution and must open all source
    resources inside the returned iterator. The iterator is closed when
    execution succeeds, fails, or is stopped early. Projection, predicates,
    and explicit row limits are applied to each emitted batch so the source
    remains correct when Polars pushes those operations into the IO plugin.

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
        source_iterator = batches(batch_size or _DEFAULT_BATCH_SIZE)
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
        is_pure=False,
    )


def register_deferred_frame_source(
    *,
    columns: Iterable[str] | None = None,
    schema: SchemaDict | None = None,
    frame: FrameFactory,
) -> pl.LazyFrame:
    """Defer a parser-backed frame factory until Polars executes the plan.

    This fallback is for resource relations whose parser already owns the
    complete domain transformation but does not expose a batch reader. Column
    names are frozen for plan construction; concrete dtypes are supplied by
    the first emitted frame. Large relations should use
    :func:`register_replayable_source` with an explicit typed schema instead.
    """

    if schema is None:
        if columns is None:
            raise ValueError("columns or schema must be provided")
        frozen_schema = dict.fromkeys(columns, pl.Unknown)
    else:
        if columns is not None:
            raise ValueError("columns and schema are mutually exclusive")
        frozen_schema = dict(schema)

    def batches(batch_size: int) -> Iterator[pl.DataFrame]:
        del batch_size
        yield frame()

    return register_replayable_source(
        schema=frozen_schema,
        batches=batches,
        validate_schema=False,
    )
