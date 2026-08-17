from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import duckdb
import polars as pl
import pytest

from bioextract.reactome import ReactomeDatabase


def _write_gmt_archive(tmp_path: Path, content: str, *, extra: bool = False) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    archive_path = tmp_path / "ReactomePathways.gmt.zip"
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr("ReactomePathways.gmt", content)
        if extra:
            archive.writestr("extra.txt", "unexpected\n")
    return archive_path


def test_gmt_parser_preserves_labels_symbols_and_human_scope(tmp_path: Path) -> None:
    archive = _write_gmt_archive(
        tmp_path,
        "GMT label with suffix\tR-HSA-1\tTP53\tgag-pol\tTP53\n"
        "Another label\tR-HSA-2\tBRCA1\n",
    )
    db = ReactomeDatabase.from_files(pathway_gene_sets=archive)
    frame = db.pathway_gene_sets().collect()
    assert frame.to_dicts() == [
        {
            "reactome_pathway_id": "R-HSA-1",
            "gene_set_name": "GMT label with suffix",
            "gene_symbol": "TP53",
        },
        {
            "reactome_pathway_id": "R-HSA-1",
            "gene_set_name": "GMT label with suffix",
            "gene_symbol": "gag-pol",
        },
        {
            "reactome_pathway_id": "R-HSA-2",
            "gene_set_name": "Another label",
            "gene_symbol": "BRCA1",
        },
    ]
    assert db.with_species("Homo sapiens").pathway_gene_sets().collect().height == 3
    assert db.with_species("Mus musculus").pathway_gene_sets().collect().columns == [
        "reactome_pathway_id",
        "gene_set_name",
        "gene_symbol",
    ]
    assert db.with_species("Mus musculus").pathway_gene_sets().collect().height == 0


def test_gmt_publication_uses_zip_provenance_and_reopens(tmp_path: Path) -> None:
    archive = _write_gmt_archive(
        tmp_path,
        "GMT label\tR-HSA-1\tTP53\n",
    )
    publication = tmp_path / "gmt.duckdb"
    source = ReactomeDatabase.from_files(
        pathway_gene_sets=archive,
        release_version="96",
    )
    result = source.write_duckdb(publication)
    assert result.tables == ("pathway_gene_set",)
    with duckdb.connect(str(publication), read_only=True) as connection:
        assert connection.execute(
            "SELECT logical_name, media_type FROM _bioextract.source_file"
        ).fetchall() == [("pathway_gene_set", "application/zip")]
        assert connection.execute(
            "SELECT value FROM _bioextract.metadata "
            "WHERE key='bioextract.source_schema_profile'"
        ).fetchone() == ("reactome-mapping-files-v5",)
    reopened = ReactomeDatabase.from_duckdb(publication)
    assert (
        reopened.pathway_gene_sets()
        .collect()
        .equals(source.pathway_gene_sets().collect())
    )


def test_gmt_archive_shape_and_label_conflicts_fail_closed(tmp_path: Path) -> None:
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="exactly one file entry"
    ):
        ReactomeDatabase.from_files(
            pathway_gene_sets=_write_gmt_archive(
                tmp_path / "extra", "label\tR-HSA-1\tTP53\n", extra=True
            )
        ).pathway_gene_sets().collect()

    conflicting = _write_gmt_archive(
        tmp_path / "conflict",
        "Label A\tR-HSA-1\tTP53\nLabel B\tR-HSA-1\tBRCA1\n",
    )
    with pytest.raises(
        (ValueError, pl.exceptions.ComputeError), match="multiple gene-set labels"
    ):
        ReactomeDatabase.from_files(
            pathway_gene_sets=conflicting
        ).pathway_gene_sets().collect()
