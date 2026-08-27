from __future__ import annotations

import tempfile
from pathlib import Path

import polars as pl
import pytest

import bioextract._publication as publication
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


def test_publication_requires_snake_case_tables(
    tmp_path: Path,
) -> None:
    dataset = _dataset(tmp_path, relation_count=2)
    with pytest.raises(ValueError, match="snake_case"):
        dataset.write_duckdb(
            tmp_path / "bad.duckdb",
            table_names={"relation_0": "relation-zero"},
        )


def test_failed_replacement_preserves_existing_publication(
    tmp_path: Path,
) -> None:
    path = tmp_path / "example.duckdb"
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
        resource_name="bad",
    )
    with pytest.raises(pl.exceptions.InvalidOperationError):
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
        resource_name="bad",
    )
    with pytest.raises(pl.exceptions.InvalidOperationError):
        failing.write_duckdb(tmp_path / "failure.duckdb")
    assert not list(tmp_path.glob("bioextract-relations-*"))


def test_unavailable_package_version_fails_before_output_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "not-created" / "publication.duckdb"

    def unavailable() -> str:
        raise RuntimeError("installed package metadata is unavailable")

    monkeypatch.setattr(publication, "require_package_version", unavailable)
    with pytest.raises(RuntimeError, match="metadata is unavailable"):
        _dataset(tmp_path).write_duckdb(destination)

    assert not destination.parent.exists()


def test_unavailable_package_version_preserves_existing_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "publication.duckdb"
    _dataset(tmp_path).write_duckdb(destination)
    original = destination.read_bytes()

    def unavailable() -> str:
        raise RuntimeError("installed package metadata is unavailable")

    monkeypatch.setattr(publication, "require_package_version", unavailable)
    with pytest.raises(RuntimeError, match="metadata is unavailable"):
        _dataset(tmp_path).write_duckdb(destination, if_exists="replace")

    assert destination.read_bytes() == original
    assert not list(tmp_path.glob(f".{destination.name}.*"))


@pytest.mark.parametrize("package_version", ["unknown", "01.0.0", "not-a-version"])
def test_package_version_must_be_resolved_canonical_pep440(
    monkeypatch: pytest.MonkeyPatch,
    package_version: str,
) -> None:
    def resolved_version(_name: str) -> str:
        return package_version

    monkeypatch.setattr(publication, "version", resolved_version)
    with pytest.raises(RuntimeError, match="package version"):
        publication.require_package_version()
