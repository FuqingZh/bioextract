from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract._tidy import TidyAsset, TidyDataset, TidySource


def test_canonical_publication_requires_final_derived_columns(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("UniProtId\nP12345\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={
            "protein_pathway": pl.DataFrame(
                {"uniprot_id": ["P12345"], "reactome_pathway_id": ["R-HSA-1"]}
            ).lazy()
        },
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="example-v1",
        source_schema_profile="example-source-v1",
        build_id_prefix="example",
        assets=(
            TidyAsset(
                "protein_pathway.parquet",
                "canonical",
                "protein_pathway",
            ),
        ),
    )

    file_duckdb = tmp_path / "example.duckdb"
    dataset.write_duckdb(file_duckdb)
    with duckdb.connect(str(file_duckdb), read_only=True) as connection:
        assert [
            row[1]
            for row in connection.execute(
                "PRAGMA table_info('protein_pathway')"
            ).fetchall()
        ] == ["uniprot_id", "reactome_pathway_id"]
        assert connection.execute("SELECT * FROM _bioextract.column_mapping").fetchall() == []


def test_derived_pascal_case_is_rejected_instead_of_normalized(tmp_path: Path) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("id\nA\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={"term": pl.DataFrame({"TermId": ["T1"]}).lazy()},
        source=TidySource("source", source, "text/plain"),
        resource_schema_version="fixture-v1",
        source_schema_profile="fixture-source-v1",
        build_id_prefix="fixture",
        assets=(TidyAsset("term.parquet", "canonical", "term"),),
    )

    with pytest.raises(ValueError, match="lowercase snake_case"):
        dataset.write_duckdb(tmp_path / "fixture.duckdb")


def test_official_headers_receive_only_required_duckdb_mapping(
    tmp_path: Path,
) -> None:
    source = tmp_path / "official.tsv"
    source.write_text("Name\tname\nA\tB\n", encoding="utf-8")
    dataset = TidyDataset(
        frames={"official": pl.DataFrame({"Name": ["A"], "name": ["B"]}).lazy()},
        source=TidySource("source", source, "text/tab-separated-values"),
        resource_schema_version="official-v1",
        source_schema_profile="official-source-v1",
        build_id_prefix="official",
        assets=(TidyAsset("official.parquet", "canonical", "official"),),
    )

    path = tmp_path / "official.duckdb"
    dataset.write_duckdb(path, source_columns={"official": ("Name", "name")})

    with duckdb.connect(str(path), read_only=True) as connection:
        assert [
            row[1]
            for row in connection.execute("PRAGMA table_info('official')").fetchall()
        ] == ["Name", "name_2"]
        assert connection.execute(
            "SELECT source_column, output_column, reason "
            "FROM _bioextract.column_mapping"
        ).fetchall() == [("name", "name_2", "case_insensitive_collision")]
