from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pytest

from bioextract.errors import CapabilityError, IntegrityError
from bioextract.kegg import KEGGDatabase


def _mapping_source(tmp_path: Path) -> KEGGDatabase:
    source = tmp_path / "hsa"
    source.mkdir(exist_ok=True)
    (source / "gene_list.tsv").write_text(
        "hsa:1\tCDS\t1..10\tGENE1; description\n", encoding="utf-8"
    )
    (source / "conv_uniprot.tsv").write_text("up:P12345\thsa:1\n", encoding="utf-8")
    (source / "gene_ko.tsv").write_text("hsa:1\tko:K00001\n", encoding="utf-8")
    (source / "gene_pathway.tsv").write_text("hsa:1\tpath:hsa00010\n", encoding="utf-8")
    return KEGGDatabase.from_mapping_files(source, organism_code="hsa")


def _brite_source(tmp_path: Path) -> KEGGDatabase:
    source = tmp_path / "brite.json"
    source.write_text(
        json.dumps(
            {
                "name": "hsa00001",
                "children": [
                    {
                        "name": "09100 Metabolism",
                        "children": [
                            {
                                "name": "09101 Carbohydrate metabolism",
                                "children": [
                                    {
                                        "name": "00010 Glycolysis [PATH:hsa00010]",
                                        "children": [
                                            {"name": "hsa:1 GENE\tK00001 enzyme"}
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return KEGGDatabase.from_brite_json(source)


def test_mapping_publication_uses_exact_v2_identity_and_three_tables(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    result = _mapping_source(tmp_path).write_duckdb(path)

    assert result.tables == ("organism", "gene_annotation", "ko_annotation")
    with KEGGDatabase.from_duckdb(path).connect() as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "2"
        assert metadata["bioextract.resource_schema_version"] == "kegg-mapping-v1.0"
        assert metadata["bioextract.source_schema_profile"] == (
            "kegg-organism-mapping-files-v2"
        )
        assert "bioextract.sources" not in metadata
        assert {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM _bioextract.table_info"
            ).fetchall()
        } == {"organism", "gene_annotation", "ko_annotation"}


def test_mapping_publication_preserves_selection_semantics(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    source = _mapping_source(tmp_path)
    expected = source.select_groups(
        {"case": ["P12345", "missing"], "repeat": ["P12345"]},
        namespace="uniprot",
    )
    source.write_duckdb(path)
    actual = KEGGDatabase.from_duckdb(path).select_groups(
        {"case": ["P12345", "missing"], "repeat": ["P12345"]},
        namespace="uniprot",
    )

    actual_matches = actual.matches().collect()
    expected_matches = expected.matches().collect()
    assert actual_matches.sort(actual_matches.columns).equals(
        expected_matches.sort(expected_matches.columns)
    )
    actual_unmatched = actual.unmatched_ids().collect()
    expected_unmatched = expected.unmatched_ids().collect()
    assert actual_unmatched.sort(actual_unmatched.columns).equals(
        expected_unmatched.sort(expected_unmatched.columns)
    )


def test_mapping_publication_rejects_table_and_capability_tampering(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    _mapping_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.table_info SET table_role='unexpected' "
            "WHERE table_name='gene_annotation'"
        )
    with pytest.raises(IntegrityError, match="row-count drift"):
        KEGGDatabase.from_duckdb(path)

    path.unlink()
    _mapping_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "DELETE FROM _bioextract.metadata WHERE key='bioextract.capability.gene_ko'"
        )
    with pytest.raises(IntegrityError, match="capability inventory"):
        KEGGDatabase.from_duckdb(path)


def test_mapping_publication_rejects_cross_species_gene(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    _mapping_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE gene_annotation SET organism_code='mmu' WHERE kegg_gene_id='hsa:1'"
        )

    with pytest.raises(IntegrityError, match="cross organism"):
        KEGGDatabase.from_duckdb(path)


def test_mapping_publication_keeps_provenance_lazy_and_atomic_destination(
    tmp_path: Path,
) -> None:
    path = tmp_path / "mapping.duckdb"
    source = _mapping_source(tmp_path)
    source.write_duckdb(path)
    before = path.read_bytes()

    with duckdb.connect(str(path), read_only=True) as connection:
        rows = connection.execute(
            "SELECT display_path, bytes, sha256 FROM _bioextract.source_file"
        ).fetchall()
    assert rows and all(
        byte_count is None and digest is None for _, byte_count, digest in rows
    )
    with pytest.raises(FileExistsError):
        source.write_duckdb(path)
    assert path.read_bytes() == before


def test_old_mapping_publication_profile_is_not_accepted(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    _mapping_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='kegg-organism-mapping-files-v1' "
            "WHERE key='bioextract.source_schema_profile'"
        )

    with pytest.raises(ValueError, match="Unsupported KEGG source schema profile"):
        KEGGDatabase.from_duckdb(path)


def test_brite_publication_contract_remains_unchanged(tmp_path: Path) -> None:
    path = tmp_path / "brite.duckdb"
    _brite_source(tmp_path).write_duckdb(path)
    reopened = KEGGDatabase.from_duckdb(path)

    with reopened.connect() as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "2"
        assert metadata["bioextract.source_schema_profile"] == "kegg-brite-json-v1"
        assert connection.execute("SELECT count(*) FROM pathway").fetchone() == (1,)
    with pytest.raises((CapabilityError, ValueError), match="BRITE"):
        reopened.select_ids(["P12345"], namespace="uniprot")


def test_brite_rejects_source_role_and_provenance_tampering(tmp_path: Path) -> None:
    path = tmp_path / "brite.duckdb"
    _brite_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.source_file SET logical_name='renamed_brite_json'"
        )
    with pytest.raises(ValueError, match="source role inventory"):
        KEGGDatabase.from_duckdb(path)

    path.unlink()
    _brite_source(tmp_path).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "INSERT INTO _bioextract.column_mapping VALUES "
            "('pathway', 'forged', 'forged', 'forged')"
        )
    with pytest.raises(ValueError, match="column provenance inventory"):
        KEGGDatabase.from_duckdb(path)
