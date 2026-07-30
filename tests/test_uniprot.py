from __future__ import annotations

import gzip
import json
import sqlite3
from pathlib import Path

import duckdb
import polars as pl
import pytest

import bioextract.uniprot._knowledgebase as knowledgebase
from bioextract.uniprot import UniProtDatabase


def _write_dat(path: Path) -> Path:
    text = """ID   TEST_HUMAN              Reviewed;         10 AA.
AC   P12345; Q11111;
DT   01-JAN-2000, sequence version 2.
DT   01-JAN-2026, entry version 9.
DE   RecName: Full=Test protein;
DE            EC=1.2.3.4;
DE   AltName: Full=First alternative;
DE   AltName: Full=Second alternative;
GN   Name=TEST;
GN   Synonyms=ALT;
GN   OrderedLocusNames=LOC1; ORFNames=ORF1;
GN   and Name=TEST2; Synonyms=ALT2;
OX   NCBI_TaxID=9606;
PE   1: Evidence at protein level;
CC   -!- FUNCTION: Demonstrates a compact fixture.
CC   -!- SUBCELLULAR LOCATION: Nucleus. Note=Fixture location.
CC   -!- ALTERNATIVE PRODUCTS:
CC       Name=Isoform 1; IsoId=P12345-1; Sequence=Displayed;
CC       Name=Isoform 2; IsoId=P12345-2; Sequence=VSP_000001;
DR   GO; GO:0003677; F:binding; EXP:UniProtKB.
DR   GeneID; 1234; -.
DR   RefSeq; NP_000001.1; NM_000001.1. [P12345-2].
DR   Ensembl; ENST1; ENSP1; ENSG1. [P12345-2].
DR   ConoServer; 1570; GIIIA [R13A].
DR   eggNOG; KOG1; Eukaryota.
KW   Reference proteome; Nucleus.
FT   VAR_SEQ         3..4
FT                   /note="AA -> GG (in isoform 2)"
FT                   /id="VSP_000001"
FT                   /evidence="ECO:0000269|PubMed:1"
FT   CONFLICT        8
FT                   /note="I -> V (in another source)"
SQ   SEQUENCE   10 AA;  1000 MW;  81FC3551E879CB1A CRC64;
     ACDEFGHIKL
//
"""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for identifier, sequence in records.items():
            handle.write(f">sp|{identifier}|fixture\n{sequence}\n")
    return path


def _write_idmapping(path: Path) -> Path:
    rows = [
        [
            "P12345",
            "TEST_HUMAN",
            "1234",
            "NP_000001.1",
            "",
            "",
            "GO:0003677",
            "",
            "",
            "",
            "",
            "",
            "9606",
            "",
            "",
            "",
            "",
            "",
            "ENSG1",
            "ENST1",
            "ENSP1",
            "",
        ]
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join("\t".join(row) for row in rows) + "\n")
    return path


def test_varsplic_identifier_lookup_uses_identifier_index() -> None:
    with sqlite3.connect(":memory:") as connection:
        knowledgebase._create_validation_index(  # pyright: ignore[reportPrivateUsage]
            connection
        )
        plan = connection.execute(
            "EXPLAIN QUERY PLAN "
            "SELECT i.primary_accession, i.isoform_id "
            "FROM isoform_identifier ii "
            "JOIN isoform i USING (primary_accession, isoform_id) "
            "WHERE ii.identifier=? AND i.sequence_status='Alternative'",
            ("P22966-1",),
        ).fetchall()

    assert any(
        "isoform_identifier_lookup" in str(detail)
        for _select_id, _order, _from, detail in plan
    )


def test_idmapping_is_lazy_and_requires_explicit_eager_scope(tmp_path: Path) -> None:
    source = _write_idmapping(tmp_path / "renamed.tab.gz")
    database = UniProtDatabase.from_idmapping(source, release_version="2026_01")

    assert database.scan_mapping().collect().height == 1
    with pytest.raises(ValueError, match="allow_all_taxa"):
        database.read_mapping()
    assert database.read_mapping(taxon_ids=["9606"]).height == 1

    result = database.write_parquet(tmp_path / "mapping.parquet", taxon_ids=["9606"])
    assert result.resource_schema_version == "uniprot-idmapping-selected-v0.1"


@pytest.mark.parametrize("compressed", [False, True])
def test_idmapping_reads_plain_and_gzip(tmp_path: Path, compressed: bool) -> None:
    path = tmp_path / ("mapping.tab.gz" if compressed else "mapping.tab")
    if compressed:
        _write_idmapping(path)
    else:
        compressed_path = _write_idmapping(tmp_path / "source.tab.gz")
        with gzip.open(compressed_path, "rt", encoding="utf-8") as source:
            path.write_text(source.read(), encoding="utf-8")
    frame = UniProtDatabase.from_idmapping(path).read_mapping(taxon_ids=["9606"])
    assert frame.select("UniProtId", "TaxId").to_dicts() == [
        {"UniProtId": "P12345", "TaxId": "9606"}
    ]
    publication = tmp_path / f"published-{compressed}.parquet"
    UniProtDatabase.from_idmapping(path).write_parquet(publication, taxon_ids=["9606"])
    with duckdb.connect() as connection:
        sources_json = connection.execute(
            "SELECT decode(value) FROM parquet_kv_metadata(?) "
            "WHERE CAST(key AS VARCHAR)='bioextract.sources'",
            [str(publication)],
        ).fetchone()
    assert sources_json is not None
    assert json.loads(sources_json[0])[0]["media_type"] == (
        "text/tab-separated-values+gzip" if compressed else "text/tab-separated-values"
    )


def test_idmapping_parquet_and_hive_inputs(tmp_path: Path) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    frame = UniProtDatabase.from_idmapping(raw).read_mapping(taxon_ids=["9606"])
    parquet = tmp_path / "mapping.parquet"
    frame.write_parquet(parquet)
    assert UniProtDatabase.from_idmapping(parquet).scan_mapping().collect().height == 1
    published = tmp_path / "republished.parquet"
    UniProtDatabase.from_idmapping(parquet).write_parquet(published, taxon_ids=["9606"])
    with duckdb.connect() as connection:
        sources_json = connection.execute(
            "SELECT decode(value) FROM parquet_kv_metadata(?) "
            "WHERE CAST(key AS VARCHAR)='bioextract.sources'",
            [str(published)],
        ).fetchone()
    assert sources_json is not None
    assert (
        json.loads(sources_json[0])[0]["media_type"] == "application/vnd.apache.parquet"
    )

    hive = tmp_path / "hive" / "TaxId=9606"
    hive.mkdir(parents=True)
    frame.write_parquet(hive / "part.parquet")
    (tmp_path / "hive" / "README.txt").write_text("ignored", encoding="utf-8")
    assert (
        UniProtDatabase.from_idmapping(tmp_path / "hive")
        .read_mapping(taxon_ids=["9606"])
        .height
        == 1
    )


def test_idmapping_schema_and_path_validation(tmp_path: Path) -> None:
    bad = tmp_path / "bad.parquet"
    pl.DataFrame({"UniProtId": ["P12345"]}).write_parquet(bad)
    with pytest.raises(ValueError, match="missing required columns"):
        UniProtDatabase.from_idmapping(bad).scan_mapping()
    with pytest.raises(FileNotFoundError):
        UniProtDatabase.from_idmapping(tmp_path / "missing.tab.gz")
    unsupported = tmp_path / "mapping.txt"
    unsupported.write_text("", encoding="utf-8")
    with pytest.raises(pl.exceptions.NoDataError, match="empty CSV"):
        UniProtDatabase.from_idmapping(unsupported).scan_mapping().collect_schema()
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no parquet"):
        UniProtDatabase.from_idmapping(empty)


def test_idmapping_publication_is_atomic(tmp_path: Path) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    destination = tmp_path / "mapping.parquet"
    destination.write_bytes(b"existing")
    database = UniProtDatabase.from_idmapping(raw)
    with pytest.raises(FileExistsError):
        database.write_parquet(destination, taxon_ids=["9606"])
    assert destination.read_bytes() == b"existing"


def test_idmapping_taxon_validation_and_all_taxa_opt_in(tmp_path: Path) -> None:
    database = UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    )
    with pytest.raises(ValueError, match="TaxId"):
        database.scan_mapping(taxon_ids=["9606", ""])
    with pytest.raises(ValueError, match="allow_all_taxa"):
        database.write_parquet(tmp_path / "blocked.parquet")
    assert database.write_parquet(
        tmp_path / "all.parquet", allow_all_taxa=True
    ).path.is_file()


def test_knowledgebase_publication_selection_and_metadata(tmp_path: Path) -> None:
    entries = _write_dat(tmp_path / "entries.bin")
    canonical = _write_fasta(tmp_path / "canonical.bin", {"P12345": "ACDEFGHIKL"})
    isoforms = _write_fasta(tmp_path / "isoforms.bin", {"P12345-2": "ACGGFGHIKL"})
    path = tmp_path / "uniprot.duckdb"

    result = UniProtDatabase.from_knowledgebase(
        entries=entries,
        canonical_sequences=canonical,
        isoform_sequences=isoforms,
        release_version="2026_01",
    ).write_duckdb(path)

    assert result.resource_schema_version == "uniprot-knowledgebase-duckdb-v1"
    database = UniProtDatabase.from_duckdb(path)
    assert database.release_version == "2026_01"
    matches = database.select_ids(["Q11111", "TEST", "missing"], namespace="uniprot")
    assert matches.extract_proteins()["UniProtId"].to_list() == ["P12345"]
    assert matches.extract_proteins()["ProteinExistence"].to_list() == [
        "1: Evidence at protein level"
    ]
    assert matches.extract_unmatched_ids()["InputId"].to_list() == [
        "TEST",
        "missing",
    ]
    assert database.select_ids(["TEST"], namespace="gene_name").extract_gene_names()[
        "GeneName"
    ].to_list() == ["TEST", "ALT", "LOC1", "ORF1", "TEST2", "ALT2"]
    assert database.select_ids(["P12345"], namespace="uniprot").extract_protein_names()[
        "ProteinName"
    ].to_list() == [
        "Test protein",
        "First alternative",
        "Second alternative",
    ]
    assert database.select_ids(["1234"], namespace="gene_id").extract_cross_references(
        databases=["GeneID"]
    )["ExternalId"].to_list() == ["1234"]
    refseq = database.select_ids(
        ["NP_000001.1"], namespace="refseq"
    ).extract_cross_references(databases=["RefSeq"])
    assert refseq.select("ExternalId", "IsoformId").to_dicts() == [
        {"ExternalId": "NP_000001.1", "IsoformId": "P12345-2"}
    ]
    con_server = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_cross_references(databases=["ConoServer"])
    assert con_server.select("ExternalId", "Properties", "IsoformId").to_dicts() == [
        {
            "ExternalId": "1570",
            "Properties": "GIIIA [R13A]",
            "IsoformId": None,
        }
    ]
    assert database.select_ids(["P12345"], namespace="uniprot").extract_sequences(
        sequence_type="all"
    ).select("SequenceType", "CRC64").to_dicts() == [
        {"SequenceType": "canonical", "CRC64": "81FC3551E879CB1A"},
        {"SequenceType": "isoform", "CRC64": None},
    ]
    selection = database.select_ids(["P12345"], namespace="uniprot")
    assert selection.extract_ec_numbers()["ECNumber"].to_list() == ["1.2.3.4"]
    assert selection.extract_go_annotations()["GOId"].to_list() == ["GO:0003677"]
    assert selection.extract_comments(comment_types=["FUNCTION"])[
        "CommentType"
    ].to_list() == ["FUNCTION"]
    assert selection.extract_subcellular_locations()[
        "SubcellularLocation"
    ].to_list() == ["Nucleus"]
    assert selection.extract_keywords()["Keyword"].to_list() == [
        "Reference proteome",
        "Nucleus",
    ]
    grouped = database.select_groups(
        {"case": ["P12345"], "control": ["missing"]},
        namespace="uniprot",
    )
    assert grouped.extract_proteins()["GroupId"].to_list() == ["case"]
    assert grouped.extract_unmatched_ids().select("GroupId", "InputId").to_dicts() == [
        {"GroupId": "control", "InputId": "missing"}
    ]
    isoforms_frame = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_isoforms()
    assert isoforms_frame.select(
        "IsoformId", "IsoformOrder", "SequenceStatus", "SequenceId"
    ).to_dicts() == [
        {
            "IsoformId": "P12345-1",
            "IsoformOrder": 1,
            "SequenceStatus": "Displayed",
            "SequenceId": "P12345:canonical",
        },
        {
            "IsoformId": "P12345-2",
            "IsoformOrder": 2,
            "SequenceStatus": "Alternative",
            "SequenceId": "P12345-2",
        },
    ]
    assert selection.extract_isoform_identifiers().select(
        "IsoformId", "Identifier", "IdentifierOrder", "IsMain"
    ).to_dicts() == [
        {
            "IsoformId": "P12345-1",
            "Identifier": "P12345-1",
            "IdentifierOrder": 1,
            "IsMain": True,
        },
        {
            "IsoformId": "P12345-2",
            "Identifier": "P12345-2",
            "IdentifierOrder": 1,
            "IsMain": True,
        },
    ]
    assert database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_isoform_variations().select(
        "IsoformId", "VariationId", "VariationOrder"
    ).to_dicts() == [
        {
            "IsoformId": "P12345-2",
            "VariationId": "VSP_000001",
            "VariationOrder": 1,
        }
    ]
    variation = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_sequence_variations()
    assert variation.select("VariationId", "Note").to_dicts() == [
        {"VariationId": "VSP_000001", "Note": "AA -> GG (in isoform 2)"}
    ]
    assert database.select_ids(
        ["P12345"], namespace="uniprot", taxon_ids=["10090"]
    ).extract_proteins().columns == [
        "GroupId",
        "InputId",
        "InputNamespace",
        "UniProtId",
        "EntryName",
        "IsReviewed",
        "TaxonId",
        "ProteinExistence",
        "SequenceLength",
        "MolecularWeight",
        "SequenceVersion",
        "EntryVersion",
    ]
    assert database.select_ids(
        [], namespace="uniprot"
    ).extract_accessions().columns == [
        "GroupId",
        "InputId",
        "InputNamespace",
        "UniProtId",
        "Accession",
        "AccessionOrder",
        "IsPrimaryAccession",
    ]
    with database.connect() as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "3"
        assert (
            metadata["bioextract.resource_schema_version"]
            == "uniprot-knowledgebase-duckdb-v1"
        )
        assert metadata["bioextract.source_schema_profile"] == "uniprotkb-flat-file-v1"
        with pytest.raises(duckdb.Error):
            connection.execute("CREATE TABLE forbidden(value INTEGER)")


def test_roles_are_validated_by_content_not_basename(tmp_path: Path) -> None:
    entries = _write_dat(tmp_path / "anything.gz")
    canonical = _write_fasta(tmp_path / "other.dat", {"P12345": "WRONG"})
    with pytest.raises(ValueError, match="Canonical FASTA"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, canonical_sequences=canonical
        ).write_duckdb(tmp_path / "bad.duckdb")

    with pytest.raises(ValueError, match="no UniProtKB records"):
        UniProtDatabase.from_knowledgebase(
            entries=_write_fasta(tmp_path / "not-dat", {"P12345": "ACD"})
        ).write_duckdb(tmp_path / "wrong-role.duckdb")


def test_release_version_is_never_inferred_from_path(tmp_path: Path) -> None:
    entries = _write_dat(tmp_path / "2026_01-uniprot.dat.gz")
    path = tmp_path / "unknown-version.duckdb"
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(path)
    assert UniProtDatabase.from_duckdb(path).release_version is None


def test_missing_required_source_field_always_fails(tmp_path: Path) -> None:
    path = tmp_path / "missing-ox.dat"
    path.write_text(
        """ID   TEST_HUMAN Reviewed; 3 AA.
AC   P12345;
SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64;
     ACD
//
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="OX"):
        UniProtDatabase.from_knowledgebase(entries=path).write_duckdb(
            tmp_path / "invalid.duckdb"
        )


def test_existing_destination_fails_before_entries_are_parsed(tmp_path: Path) -> None:
    entries = tmp_path / "invalid.dat"
    entries.write_text("not a UniProt record\n", encoding="utf-8")
    destination = tmp_path / "existing.duckdb"
    destination.write_bytes(b"existing")

    with pytest.raises(FileExistsError):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(destination)
    assert destination.read_bytes() == b"existing"


def test_accessions_must_be_unique_across_records(tmp_path: Path) -> None:
    entries = tmp_path / "duplicate-accession.dat"
    entries.write_text(
        """ID   FIRST_HUMAN Reviewed; 3 AA.
AC   P11111; Q99999;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64;
     ACD
//
ID   SECOND_HUMAN Reviewed; 3 AA.
AC   P22222; Q99999;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  300 MW;  69CB1DB000000000 CRC64;
     AEF
//
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="reused across records: Q99999"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "duplicate.duckdb"
        )


def test_metadata_v3_requires_validation_issue_table(tmp_path: Path) -> None:
    path = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute("DROP TABLE _bioextract.validation_issue")
    with pytest.raises(ValueError, match="validation_issue"):
        UniProtDatabase.from_duckdb(path)


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        (
            "ALTER TABLE protein ALTER sequence_length TYPE VARCHAR",
            "physical schema mismatch",
        ),
        (
            "ALTER TABLE protein RENAME sequence_length TO reported_length",
            "physical schema mismatch",
        ),
        ("CREATE TABLE unexpected(value INTEGER)", "relation inventory mismatch"),
        (
            "INSERT INTO _bioextract.table_info VALUES ('unexpected', 'canonical', 0)",
            "table_info inventory mismatch",
        ),
        (
            "UPDATE _bioextract.metadata SET value='unsupported-profile' "
            "WHERE key='bioextract.source_schema_profile'",
            "source schema profile",
        ),
    ],
)
def test_from_duckdb_rejects_physical_contract_corruption(
    tmp_path: Path, corruption: str, message: str
) -> None:
    path = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(path)
    with duckdb.connect(str(path)) as connection:
        connection.execute(corruption)

    with pytest.raises(ValueError, match=message):
        UniProtDatabase.from_duckdb(path)


def test_selection_and_fasta_inputs_reject_empty_values(tmp_path: Path) -> None:
    path = tmp_path / "uniprot.duckdb"
    entries = _write_dat(tmp_path / "entries.dat.gz")
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(path)
    database = UniProtDatabase.from_duckdb(path)
    with pytest.raises(ValueError, match="TaxId"):
        database.select_ids(["P12345"], namespace="uniprot", taxon_ids=[""])
    with pytest.raises(ValueError, match="group labels"):
        database.select_groups({"": ["P12345"]}, namespace="uniprot")
    with pytest.raises(ValueError, match="unique after normalization"):
        database.select_groups(
            {" A ": ["P12345"], "A": ["Q11111"]}, namespace="uniprot"
        )

    empty_fasta = tmp_path / "empty-sequence.fasta"
    empty_fasta.write_text(">sp|P12345|TEST_HUMAN\n   \n\t\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid sequence characters"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, canonical_sequences=empty_fasta
        ).write_duckdb(tmp_path / "invalid-fasta.duckdb")

    embedded_blank = tmp_path / "embedded-blank.fasta"
    embedded_blank.write_text(">sp|P12345-2|TEST_HUMAN\nACD\n\nEF\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid sequence characters"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, isoform_sequences=embedded_blank
        ).write_duckdb(tmp_path / "embedded-blank.duckdb")

    invalid_fasta = tmp_path / "invalid-sequence.fasta"
    invalid_fasta.write_text(">sp|P12345-2|TEST_HUMAN\nAC D1!\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid sequence characters"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, isoform_sequences=invalid_fasta
        ).write_duckdb(tmp_path / "invalid-characters.duckdb")


def test_dat_crc64_must_match_sequence(tmp_path: Path) -> None:
    entries = tmp_path / "bad-checksum.dat"
    entries.write_text(
        """ID   TEST_HUMAN Reviewed; 3 AA.
AC   P12345;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  300 MW;  0000000000000000 CRC64;
     ACD
//
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="CRC64 mismatch"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "bad-checksum.duckdb"
        )


def test_cross_entry_external_isoform_contexts_and_owner_materialization(
    tmp_path: Path,
) -> None:
    entries = tmp_path / "paired-entries.dat"
    entries.write_text(
        """ID   7UP1_DROME Reviewed; 3 AA.
AC   P16375;
OX   NCBI_TaxID=7227;
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative splicing; Named isoforms=3;
CC       Name=B; IsoId=P16375-1; Sequence=Displayed;
CC       Name=C; IsoId=P16375-2, P22966-1; Sequence=VSP_013348;
CC       Name=A; IsoId=P16376-1; Sequence=External;
FT   VAR_SEQ         2
FT                   /note="C -> G (in isoform C)"
FT                   /id="VSP_013348"
SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64;
     ACD
//
ID   7UP2_DROME Reviewed; 3 AA.
AC   P16376;
OX   NCBI_TaxID=7227;
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative splicing; Named isoforms=3;
CC       Name=A; IsoId=P16376-1; Sequence=Displayed;
CC       Name=B; IsoId=P16375-1; Sequence=External;
CC       Name=C; IsoId=P16375-2, P22966-1; Sequence=External;
SQ   SEQUENCE   3 AA;  300 MW;  69CB1DB000000000 CRC64;
     AEF
//
""",
        encoding="utf-8",
    )
    path = tmp_path / "paired.duckdb"
    varsplic = _write_fasta(tmp_path / "paired-varsplic.fasta.gz", {"P22966-1": "AGD"})
    UniProtDatabase.from_knowledgebase(
        entries=entries, isoform_sequences=varsplic
    ).write_duckdb(path)

    with UniProtDatabase.from_duckdb(path).connect() as connection:
        assert connection.execute(
            "SELECT primary_accession, isoform_id, sequence_status, sequence_id "
            "FROM protein_isoform ORDER BY primary_accession, isoform_id"
        ).fetchall() == [
            ("P16375", "P16375-1", "Displayed", "P16375:canonical"),
            ("P16375", "P16375-2", "Alternative", "P16375-2"),
            ("P16375", "P16376-1", "External", None),
            ("P16376", "P16375-1", "External", None),
            ("P16376", "P16375-2", "External", None),
            ("P16376", "P16376-1", "Displayed", "P16376:canonical"),
        ]
        assert connection.execute(
            "SELECT primary_accession, isoform_id, identifier, identifier_order, "
            "is_main FROM protein_isoform_identifier "
            "WHERE isoform_id='P16375-2' "
            "ORDER BY primary_accession, identifier_order"
        ).fetchall() == [
            ("P16375", "P16375-2", "P16375-2", 1, True),
            ("P16375", "P16375-2", "P22966-1", 2, False),
            ("P16376", "P16375-2", "P16375-2", 1, True),
            ("P16376", "P16375-2", "P22966-1", 2, False),
        ]
        assert connection.execute(
            "SELECT primary_accession, isoform_id, variation_id "
            "FROM protein_isoform_variation"
        ).fetchall() == [("P16375", "P16375-2", "VSP_013348")]
        assert connection.execute(
            "SELECT sequence_id, primary_accession, sequence_type, sequence "
            "FROM protein_sequence WHERE sequence_type='isoform'"
        ).fetchall() == [("P16375-2", "P16375", "isoform", "AGD")]
        assert connection.execute(
            "SELECT ii.identifier, i.sequence_id "
            "FROM protein_isoform_identifier ii "
            "JOIN protein_isoform i USING (primary_accession, isoform_id) "
            "WHERE ii.primary_accession='P16375' "
            "AND ii.isoform_id='P16375-2' ORDER BY ii.identifier_order"
        ).fetchall() == [
            ("P16375-2", "P16375-2"),
            ("P22966-1", "P16375-2"),
        ]

    database = UniProtDatabase.from_duckdb(path)
    assert database.select_ids(["P22966-1"], namespace="isoform_id").extract_proteins()[
        "UniProtId"
    ].to_list() == ["P16375", "P16376"]
    assert database.select_ids(
        ["P16375"], namespace="uniprot"
    ).extract_isoform_variations()["VariationId"].to_list() == ["VSP_013348"]
    assert (
        database.select_ids(["P16376"], namespace="uniprot")
        .extract_isoform_variations()
        .is_empty()
    )
