from __future__ import annotations

from pathlib import Path

import duckdb
import pytest

from bioextract.errors import IntegrityError
from bioextract.reactome import ReactomeDatabase


def _write_all_levels_source(tmp_path: Path) -> Path:
    source = tmp_path / "UniProt2Reactome_All_Levels.txt"
    source.write_text(
        "P04637\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway\tTAS\tHomo sapiens\n",
        encoding="utf-8",
    )
    return source


def test_all_level_publication_has_exact_role_inventory(tmp_path: Path) -> None:
    source = ReactomeDatabase.from_files(
        uniprot_all_levels=_write_all_levels_source(tmp_path),
        release_version="96",
    )
    publication = tmp_path / "reactome.duckdb"
    source.write_duckdb(publication)

    with duckdb.connect(str(publication), read_only=True) as connection:
        roles = connection.execute(
            "SELECT logical_name FROM _bioextract.source_file"
        ).fetchall()
        assert roles == [("uniprot_pathway_all_level",)]
        tables = connection.execute(
            "SELECT table_name FROM _bioextract.table_info"
        ).fetchall()
        assert tables == [("uniprot_pathway_all_level",)]

    assert ReactomeDatabase.from_duckdb(publication).release_version == "96"


def test_reopen_rejects_pre_v02_reactome_profile(tmp_path: Path) -> None:
    publication = tmp_path / "reactome.duckdb"
    ReactomeDatabase.from_files(
        uniprot_all_levels=_write_all_levels_source(tmp_path)
    ).write_duckdb(publication)
    with duckdb.connect(str(publication)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata "
            "SET value='reactome-mapping-files-v1' "
            "WHERE key='bioextract.source_schema_profile'"
        )

    with pytest.raises(IntegrityError, match="source schema profile"):
        ReactomeDatabase.from_duckdb(publication)
