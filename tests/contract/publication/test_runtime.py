from pathlib import Path
from typing import Any

import duckdb
import polars as pl
import pytest

from bioextract import _publication
from bioextract._tidy import TidyAsset, TidyDataset, TidySource


def test_publication_connections_share_the_polars_thread_budget(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={"relation": pl.DataFrame({"id": ["T1"]}).lazy()},
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=(TidyAsset("relation.parquet", "canonical", "relation"),),
        resource_name="example",
    )
    original_connect = duckdb.connect
    calls: list[tuple[bool, dict[str, Any] | None]] = []

    def connect(
        database: str | Path = ":memory:",
        read_only: bool = False,
        config: dict[str, Any] | None = None,
    ) -> duckdb.DuckDBPyConnection:
        calls.append((read_only, config))
        return original_connect(database=database, read_only=read_only, config=config)

    monkeypatch.setattr(_publication.duckdb, "connect", connect)

    dataset.write_duckdb(tmp_path / "publication.duckdb")

    expected_config = {"threads": str(pl.thread_pool_size())}
    assert calls == [(False, expected_config), (True, expected_config)]
