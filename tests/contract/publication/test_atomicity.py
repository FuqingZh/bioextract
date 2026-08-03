from __future__ import annotations

import tempfile
from pathlib import Path

import polars as pl
import pytest

from bioextract._tidy import TidyAsset, TidyDataset, TidySource


def _dataset(tmp_path: Path, *, relation_count: int = 1) -> TidyDataset:
    source = tmp_path / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    frames = {
        f"relation_{index}": pl.DataFrame({"id": [f"T{index}"]}).lazy()
        for index in range(relation_count)
    }
    return TidyDataset(
        frames=frames,
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=tuple(
            TidyAsset(
                path=f"relation_{index}.parquet",
                kind="canonical",
                frame_name=f"relation_{index}",
            )
            for index in range(relation_count)
        ),
        resource_name="example",
        release_version="2026-07-29",
    )


def test_publication_is_atomic_and_requires_snake_case_tables(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.parquet"
    _dataset(tmp_path).write_parquet(path)
    original = path.read_bytes()

    with pytest.raises(FileExistsError):
        _dataset(tmp_path).write_parquet(path)
    assert path.read_bytes() == original

    dataset = _dataset(tmp_path, relation_count=2)
    with pytest.raises(ValueError, match="snake_case"):
        dataset.write_duckdb(
            tmp_path / "bad.duckdb",
            table_names={"relation_0": "relation-zero"},
        )


@pytest.mark.parametrize("container", ["parquet", "duckdb"])
def test_failed_replacement_preserves_existing_publication(
    tmp_path: Path,
    container: str,
) -> None:
    path = tmp_path / f"example.{container}"
    if container == "parquet":
        _dataset(tmp_path).write_parquet(path)
    else:
        _dataset(tmp_path).write_duckdb(path)
    original = path.read_bytes()

    source = tmp_path / "bad-source.tsv"
    source.write_text("value\nbad\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={
            "relation": pl.DataFrame({"value": ["bad"]})
            .lazy()
            .select(pl.col("value").cast(pl.Int64))
        },
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="bad-v1",
        source_schema_profile="bad-source-v1",
        build_id_prefix="bad",
        assets=(TidyAsset("relation.parquet", "canonical", "relation"),),
    )
    with pytest.raises(pl.exceptions.InvalidOperationError):
        if container == "parquet":
            dataset.write_parquet(path, if_exists="replace")
        else:
            dataset.write_duckdb(path, if_exists="replace")

    assert path.read_bytes() == original


def test_duckdb_transfer_parquet_directory_is_cleaned_on_success_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    _dataset(tmp_path).write_duckdb(tmp_path / "success.duckdb")
    assert not list(tmp_path.glob("bioextract-relations-*"))

    source = tmp_path / "bad-transfer-source.tsv"
    source.write_text("value\nbad\n", encoding="utf-8")
    failing = TidyDataset(
        frames={
            "relation": pl.DataFrame({"value": ["bad"]})
            .lazy()
            .select(pl.col("value").cast(pl.Int64))
        },
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="bad-v1",
        source_schema_profile="bad-source-v1",
        build_id_prefix="bad",
        assets=(TidyAsset("relation.parquet", "canonical", "relation"),),
    )
    with pytest.raises(pl.exceptions.InvalidOperationError):
        failing.write_duckdb(tmp_path / "failure.duckdb")
    assert not list(tmp_path.glob("bioextract-relations-*"))
