from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import polars as pl

from bioextract._lazy import register_replayable_source
from bioextract._publication import (
    RelationSpec,
    SourceFileRecord,
    write_duckdb_publication,
)


class _RequestLike(Protocol):
    @property
    def columns(self) -> tuple[str, ...] | None: ...

    @property
    def predicate(self) -> pl.Expr | None: ...

    @property
    def n_rows(self) -> int | None: ...

    @property
    def batch_size(self) -> int | None: ...

    @property
    def effective_batch_size(self) -> int: ...


@dataclass
class _Probe:
    batches: tuple[pl.DataFrame, ...]
    requests: list[_RequestLike]
    opened: int = 0
    closed: int = 0

    def source(self, request: _RequestLike) -> Iterator[pl.DataFrame]:
        self.requests.append(request)
        self.opened += 1
        try:
            for frame in self.batches:
                columns = request.columns
                yield frame if columns is None else frame.select(columns)
        finally:
            self.closed += 1


def _probe() -> _Probe:
    return _Probe(
        requests=[],
        batches=(
            pl.DataFrame({"id": [1, 2], "kind": ["keep", "drop"], "value": [10, 20]}),
            pl.DataFrame({"id": [3, 4], "kind": ["keep", "keep"], "value": [30, 40]}),
        ),
    )


def _frame(probe: _Probe, *, is_pure: bool = False) -> pl.LazyFrame:
    return register_replayable_source(
        schema={"id": pl.Int64, "kind": pl.String, "value": pl.Int64},
        batches=probe.source,
        is_pure=is_pure,
    )


def test_request_contains_projection_and_predicate_columns() -> None:
    probe = _probe()

    result = _frame(probe).filter(pl.col("kind") == "keep").select("id").collect()

    assert result.to_dicts() == [{"id": 1}, {"id": 3}, {"id": 4}]
    assert probe.opened == 1
    assert probe.closed == 1
    assert len(probe.requests) == 1
    assert probe.requests[0].columns == ("id", "kind")


def test_projection_is_sent_to_the_producer() -> None:
    probe = _probe()

    result = _frame(probe).select("value").collect()

    assert result.to_dicts() == [
        {"value": 10},
        {"value": 20},
        {"value": 30},
        {"value": 40},
    ]
    assert probe.requests[0].columns == ("value",)


def test_head_stops_after_the_requested_rows_and_closes_source() -> None:
    probe = _probe()

    result = _frame(probe).head(1).collect()

    assert result.to_dicts() == [{"id": 1, "kind": "keep", "value": 10}]
    assert probe.requests[0].n_rows == 1
    assert probe.closed == 1


def test_zero_row_slice_does_not_open_the_producer() -> None:
    probe = _probe()

    assert _frame(probe).head(0).collect().is_empty()

    assert probe.opened == 0
    assert probe.closed == 0


def test_repeated_collection_reopens_a_replayable_source() -> None:
    probe = _probe()
    frame = _frame(probe)

    assert frame.collect().height == 4
    assert frame.collect().height == 4

    assert probe.opened == 2
    assert probe.closed == 2


def test_pure_sources_are_shared_by_collect_all() -> None:
    probe = _probe()
    frame = _frame(probe, is_pure=True)

    first, second = pl.collect_all([frame.select("id"), frame.select("value")])

    assert first.height == 4
    assert second.height == 4
    assert probe.opened == 1
    assert probe.closed == 1


def test_publication_writer_consumes_one_shot_batches_without_prebuffer(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "source.txt"
    source_path.write_text("fixture\n", encoding="utf-8")
    emitted_sizes: list[int] = []

    def batches(_request: _RequestLike) -> Iterator[pl.DataFrame]:
        for values in ([1, 2], [3, 4]):
            emitted_sizes.append(len(values))
            yield pl.DataFrame({"id": values})

    frame = register_replayable_source(
        schema={"id": pl.Int64},
        batches=batches,
    )
    output = tmp_path / "publication.duckdb"
    result = write_duckdb_publication(
        [RelationSpec("items", frame)],
        output,
        resource_name="fixture",
        resource_schema_version="fixture-v1",
        source_schema_profile="fixture-v1",
        sources=(
            SourceFileRecord(
                "source",
                source_path,
                "text/plain",
                source_path.stat().st_size,
            ),
        ),
    )

    assert result.row_counts == {"items": 4}
    assert emitted_sizes == [2, 2]
