from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from bioextract.errors import CapabilityError, IntegrityError
from bioextract.kegg import KEGGDatabase


def write_mapping_tree(tmp_path: Path) -> Path:
    root = tmp_path / "mapping"
    for code, accession, ncbi, ko in (
        ("hsa", "P12345", "101", "K00001"),
        ("mmu", "P12345", "201", "K00002"),
    ):
        organism = root / code
        organism.mkdir(parents=True)
        (organism / "gene_list.tsv").write_text(
            f"{code}:1\tCDS\t1..100\tGENE1, ALIAS2, ALIAS1; alpha description\n"
            f"{code}:2\tCDS\tcomplement(200..300)\tsecond description\n",
            encoding="utf-8",
        )
        (organism / "conv_uniprot.tsv").write_text(
            f"up:{accession}\t{code}:1\nup:Q9Y243\t{code}:2\n",
            encoding="utf-8",
        )
        (organism / "conv_ncbi_geneid.tsv").write_text(
            f"ncbi-geneid:{ncbi}\t{code}:1\n",
            encoding="utf-8",
        )
        (organism / "gene_ko.tsv").write_text(
            f"{code}:1\tko:{ko}\n{code}:2\tko:K00003\n",
            encoding="utf-8",
        )
        (organism / "gene_pathway.tsv").write_text(
            f"{code}:1\tpath:{code}00010\n",
            encoding="utf-8",
        )
    (root / "organism").mkdir()
    (root / "organism" / "list_organism.tsv").write_text(
        "T01001\thsa\tHomo sapiens\tEukaryotes;Animals\n"
        "T01002\tmmu\tMus musculus\tEukaryotes;Animals\n",
        encoding="utf-8",
    )
    (root / "ko").mkdir()
    (root / "ko" / "ko_pathway.tsv").write_text(
        "ko:K00001\tpath:map00010\nko:K00001\tpath:ko00010\nko:K00002\tpath:map00020\n",
        encoding="utf-8",
    )
    return root


def test_directory_factory_exposes_fixed_nested_relations_lazily(
    tmp_path: Path,
) -> None:
    db = KEGGDatabase.from_mapping_directory(write_mapping_tree(tmp_path))

    gene_lf = db.gene_annotations()
    assert isinstance(gene_lf, pl.LazyFrame)
    assert gene_lf.collect_schema()["ko_mappings"] == pl.List(
        pl.Struct({"ko_id": pl.String})
    )
    genes = gene_lf.collect().sort("organism_code", "kegg_gene_id")
    assert genes.height == 4
    first = genes.row(0, named=True)
    assert first["gene_aliases"] == ["ALIAS1", "ALIAS2"]
    assert first["uniprot_mappings"] == [{"uniprot_id": "P12345"}]
    assert first["pathway_mappings"] == [
        {"kegg_pathway_id": "hsa00010", "pathway_map_id": "map00010"}
    ]
    assert genes.row(1, named=True)["gene_symbol"] is None
    assert genes.row(1, named=True)["gene_description"] == "second description"

    organisms = db.organisms().collect().sort("organism_code")
    assert organisms.select("organism_code", "genome_id").to_dicts() == [
        {"organism_code": "hsa", "genome_id": "T01001"},
        {"organism_code": "mmu", "genome_id": "T01002"},
    ]
    assert organisms["taxonomy_lineage"].to_list() == [
        ["Eukaryotes", "Animals"],
        ["Eukaryotes", "Animals"],
    ]


def test_scope_prunes_organisms_but_keeps_global_kos(tmp_path: Path) -> None:
    db = KEGGDatabase.from_mapping_directory(write_mapping_tree(tmp_path))
    hsa = db.with_organisms(["hsa", "hsa"])

    assert hsa.organisms().collect()["organism_code"].to_list() == ["hsa"]
    assert hsa.gene_annotations().collect()["organism_code"].unique().to_list() == [
        "hsa"
    ]
    assert hsa.ko_annotations().collect().sort("ko_id")["ko_id"].to_list() == [
        "K00001",
        "K00002",
        "K00003",
    ]
    with pytest.raises(ValueError, match="without case rewriting"):
        db.with_organisms(["HSA"])


def test_direct_and_via_ko_pathways_remain_distinct(tmp_path: Path) -> None:
    hsa = KEGGDatabase.from_mapping_directory(
        write_mapping_tree(tmp_path)
    ).with_organisms(["hsa"])

    direct = hsa.gene_pathways().collect().filter(pl.col("kegg_gene_id") == "hsa:1")
    assert direct["pathway_mappings"].to_list() == [
        [{"kegg_pathway_id": "hsa00010", "pathway_map_id": "map00010"}]
    ]
    via = hsa.gene_pathways_via_ko().collect()
    assert via.row(0, named=True)["pathway_mappings"] == [
        {
            "ko_id": "K00001",
            "kegg_pathway_id": "ko00010",
            "pathway_namespace": "ko",
            "pathway_map_id": "map00010",
        },
        {
            "ko_id": "K00001",
            "kegg_pathway_id": "map00010",
            "pathway_namespace": "map",
            "pathway_map_id": "map00010",
        },
    ]


def test_selection_retains_many_species_and_nested_input_lineage(
    tmp_path: Path,
) -> None:
    db = KEGGDatabase.from_mapping_directory(write_mapping_tree(tmp_path))
    selection = db.select_groups(
        {
            "case": ["sp|P12345|EXAMPLE", "up:P12345", "missing"],
            "control": ["Q9Y243"],
        },
        namespace="uniprot",
    )

    matches = selection.matches().collect().sort("group_id", "organism_code")
    assert matches.filter(pl.col("input_id") == "P12345").select(
        "organism_code", "kegg_gene_id"
    ).to_dicts() == [
        {"organism_code": "hsa", "kegg_gene_id": "hsa:1"},
        {"organism_code": "mmu", "kegg_gene_id": "mmu:1"},
    ]
    annotations = selection.gene_annotations().collect()
    assert "inputs" in annotations.columns
    assert annotations.filter(pl.col("kegg_gene_id") == "hsa:1").row(0, named=True)[
        "inputs"
    ] == [
        {
            "group_id": "case",
            "input_id": "P12345",
            "input_namespace": "uniprot",
        }
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "case", "input_id": "missing"}
    ]


def test_partial_file_factory_uses_null_and_empty_capability_semantics(
    tmp_path: Path,
) -> None:
    source = tmp_path / "hsa"
    source.mkdir()
    (source / "gene_list.tsv").write_text(
        "hsa:1\tCDS\t1\tGENE1; description\n", encoding="utf-8"
    )
    (source / "gene_ko.tsv").write_text("", encoding="utf-8")
    db = KEGGDatabase.from_mapping_files(source, organism_code="hsa")

    row = db.gene_annotations().collect().row(0, named=True)
    assert row["gene_aliases"] == []
    assert row["ko_mappings"] == []
    assert row["uniprot_mappings"] is None
    with pytest.raises(CapabilityError, match="gene_pathway"):
        db.gene_pathways()
    with pytest.raises(CapabilityError, match="uniprot_conversion"):
        db.select_ids(["P12345"], namespace="uniprot").matches()


def test_strict_parser_reports_role_path_and_line(tmp_path: Path) -> None:
    root = write_mapping_tree(tmp_path)
    (root / "hsa" / "gene_ko.tsv").write_text("mmu:1\tko:K00001\n", encoding="utf-8")
    db = KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"])

    with pytest.raises(
        (IntegrityError, pl.exceptions.ComputeError),
        match=r"gene_ko.*line=1.*does not match",
    ):
        db.gene_annotations().collect()


def test_publication_has_three_tables_capabilities_and_source_parity(
    tmp_path: Path,
) -> None:
    source = KEGGDatabase.from_mapping_directory(
        write_mapping_tree(tmp_path), release_version="2026-06"
    ).with_organisms(["hsa"])
    path = tmp_path / "mapping.duckdb"
    result = source.write_duckdb(path)

    assert result.tables == ("organism", "gene_annotation", "ko_annotation")
    reopened = KEGGDatabase.from_duckdb(path)
    assert (
        reopened.gene_annotations()
        .collect()
        .equals(source.gene_annotations().collect())
    )
    assert (
        reopened.ko_annotations()
        .collect()
        .sort("ko_id")
        .equals(source.ko_annotations().collect().sort("ko_id"))
    )
    with reopened.connect() as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.resource_schema_version"] == "kegg-mapping-v1.0"
        assert metadata["bioextract.source_schema_profile"] == (
            "kegg-organism-mapping-files-v2"
        )
        assert metadata["bioextract.organism_scope_mode"] == "selected"
        assert metadata["bioextract.release_version"] == "2026-06"
        assert metadata["bioextract.capability.ko_pathway"] == "true"
        sources = connection.execute(
            "SELECT logical_name, display_path, bytes, sha256 "
            "FROM _bioextract.source_file ORDER BY logical_name"
        ).fetchall()
        assert len(sources) == 7
        assert all(
            byte_count is None and digest is None
            for _, _, byte_count, digest in sources
        )


def test_lazy_plan_construction_does_not_enumerate_or_open_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = write_mapping_tree(tmp_path)
    database = KEGGDatabase.from_mapping_directory(root)

    def forbidden(*args: object, **kwargs: object) -> object:
        del args, kwargs
        raise AssertionError("lazy plan construction touched source data")

    monkeypatch.setattr(Path, "iterdir", forbidden)
    monkeypatch.setattr(Path, "read_bytes", forbidden)
    frame = database.gene_annotations()
    assert frame.collect_schema()["kegg_gene_id"] == pl.String


def test_mapping_publisher_does_not_read_sources_for_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = write_mapping_tree(tmp_path)
    database = KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"])

    def forbidden(*args: object, **kwargs: object) -> object:
        del args, kwargs
        raise AssertionError("provenance unexpectedly materialized source bytes")

    monkeypatch.setattr(Path, "read_bytes", forbidden)
    path = tmp_path / "once.duckdb"
    database.write_duckdb(path)
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT count(*) FROM _bioextract.source_file "
            "WHERE bytes IS NULL AND sha256 IS NULL"
        ).fetchone() == (7,)


def test_mapping_publication_failure_cleans_stage_and_preserves_target(
    tmp_path: Path,
) -> None:
    root = write_mapping_tree(tmp_path)
    valid = KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"])
    target = tmp_path / "mapping.duckdb"
    valid.write_duckdb(target)
    before = target.read_bytes()
    (root / "hsa" / "gene_ko.tsv").write_text("mmu:1\tko:K00001\n", encoding="utf-8")

    with pytest.raises(IntegrityError, match="cause=invalid_required_field"):
        valid.write_duckdb(target, if_exists="replace")

    assert target.read_bytes() == before
    assert not list(tmp_path.glob(".mapping.duckdb.*.duckdb"))
    assert not list(tmp_path.glob("bioextract-kegg-duckdb-*"))


def test_native_writer_preserves_zero_byte_roles_and_deduplicates_edges(
    tmp_path: Path,
) -> None:
    root = write_mapping_tree(tmp_path)
    gene_ko = root / "hsa" / "gene_ko.tsv"
    gene_ko.write_text("hsa:1\tko:K00001\nhsa:1\tko:K00001\n", encoding="utf-8")
    (root / "hsa" / "conv_uniprot.tsv").write_text("", encoding="utf-8")
    path = tmp_path / "mapping.duckdb"

    KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"]).write_duckdb(path)

    row = (
        KEGGDatabase.from_duckdb(path)
        .gene_annotations()
        .collect()
        .filter(pl.col("kegg_gene_id") == "hsa:1")
        .row(0, named=True)
    )
    assert row["ko_mappings"] == [{"ko_id": "K00001"}]
    with duckdb.connect(str(path), read_only=True) as connection:
        assert connection.execute(
            "SELECT bytes, sha256 FROM _bioextract.source_file "
            "WHERE logical_name='organism/hsa/uniprot_conversion'"
        ).fetchone() == (None, None)
        assert not connection.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_name LIKE '_kegg_%'"
        ).fetchall()


@pytest.mark.parametrize(
    "payload",
    [b"hsa:1\n", b"hsa:1\tko:K00001\ttoo-many\n"],
)
def test_native_writer_rejects_malformed_role_rows(
    tmp_path: Path,
    payload: bytes,
) -> None:
    root = write_mapping_tree(tmp_path)
    (root / "hsa" / "gene_ko.tsv").write_bytes(payload)

    with pytest.raises(
        IntegrityError,
        match=r"role='gene_ko'.*cause=(engine_parse|csv_rejects)",
    ):
        KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"]).write_duckdb(
            tmp_path / "mapping.duckdb"
        )


def test_native_writer_rejects_conflicting_metadata_after_duplicate_dedup(
    tmp_path: Path,
) -> None:
    root = write_mapping_tree(tmp_path)
    (root / "hsa" / "gene_list.tsv").write_text(
        "hsa:1\tCDS\t1..100\tGENE1; first\nhsa:1\tCDS\t1..100\tGENE1; second\n",
        encoding="utf-8",
    )

    with pytest.raises(IntegrityError, match="role='gene_list'.*conflicting_metadata"):
        KEGGDatabase.from_mapping_directory(root).with_organisms(["hsa"]).write_duckdb(
            tmp_path / "mapping.duckdb"
        )


def test_reopened_publication_can_narrow_but_not_expand(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    KEGGDatabase.from_mapping_directory(write_mapping_tree(tmp_path)).write_duckdb(path)
    reopened = KEGGDatabase.from_duckdb(path)

    assert (
        reopened.with_organisms(["mmu"])
        .organisms()
        .collect()
        .to_dicts()[0]["organism_code"]
        == "mmu"
    )
    with pytest.raises(CapabilityError, match="does not contain"):
        reopened.with_organisms(["eco"])


def test_explicit_role_does_not_fall_back_and_duplicate_paths_fail(
    tmp_path: Path,
) -> None:
    source = tmp_path / "hsa"
    source.mkdir()
    role = source / "gene_ko.tsv"
    role.write_text("hsa:1\tko:K00001\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        KEGGDatabase.from_mapping_files(
            source,
            organism_code="hsa",
            gene_ko=source / "missing.tsv",
        )
    with pytest.raises(ValueError, match="multiple roles"):
        KEGGDatabase.from_mapping_files(
            organism_code="hsa",
            gene_ko=role,
            gene_pathway=role,
        )


def test_publication_rejects_capability_tampering(tmp_path: Path) -> None:
    path = tmp_path / "mapping.duckdb"
    KEGGDatabase.from_mapping_directory(write_mapping_tree(tmp_path)).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            "UPDATE _bioextract.metadata SET value='maybe' "
            "WHERE key='bioextract.capability.gene_ko'"
        )

    with pytest.raises(IntegrityError, match="capabilities must be true or false"):
        KEGGDatabase.from_duckdb(path)
