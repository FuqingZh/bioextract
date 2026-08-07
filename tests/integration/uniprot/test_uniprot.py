from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import polars as pl
import pytest

import bioextract.uniprot._query as uniprot_query
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
CC   -!- FUNCTION: Demonstrates a mem-
CC       brane fixture.
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
SQ   SEQUENCE   10 AA;  1132 MW;  81FC3551E879CB1A CRC64;
     ACDEFGHIKL
//
"""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


def _write_fasta(
    path: Path, records: dict[str, str], *, compressed: bool = True
) -> Path:
    opener = gzip.open if compressed else Path.open
    with opener(path, "wt", encoding="utf-8") as handle:
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


def test_idmapping_is_lazy_and_requires_explicit_eager_scope(tmp_path: Path) -> None:
    source = _write_idmapping(tmp_path / "renamed.tab.gz")
    database = UniProtDatabase.from_idmapping(source, release_version="2026_01")

    assert database.scan_mapping().collect().height == 1
    with pytest.raises(ValueError, match="allow_all_taxa"):
        database.read_mapping()
    assert database.read_mapping(taxon_ids=["9606"]).height == 1

    result = database.write_duckdb(tmp_path / "mapping.duckdb", taxon_ids=["9606"])
    assert result.resource_schema_version == "uniprot-idmapping-duckdb-v1"


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
    assert frame.select("uniprot_id", "tax_id").to_dicts() == [
        {"uniprot_id": "P12345", "tax_id": "9606"}
    ]
    publication = tmp_path / f"published-{compressed}.duckdb"
    UniProtDatabase.from_idmapping(path).write_duckdb(publication, taxon_ids=["9606"])
    with UniProtDatabase.from_duckdb(publication).connect() as connection:
        source_media_type = connection.execute(
            "SELECT media_type FROM _bioextract.source_file"
        ).fetchone()
    assert source_media_type == (
        "text/tab-separated-values+gzip" if compressed else "text/tab-separated-values",
    )


def test_idmapping_parquet_and_hive_inputs(tmp_path: Path) -> None:
    raw = _write_idmapping(tmp_path / "mapping.tab.gz")
    frame = UniProtDatabase.from_idmapping(raw).read_mapping(taxon_ids=["9606"])
    parquet = tmp_path / "mapping.parquet"
    frame.write_parquet(parquet)
    assert UniProtDatabase.from_idmapping(parquet).scan_mapping().collect().height == 1
    published = tmp_path / "republished.duckdb"
    UniProtDatabase.from_idmapping(parquet).write_duckdb(published, taxon_ids=["9606"])
    with UniProtDatabase.from_duckdb(published).connect() as connection:
        assert connection.execute(
            "SELECT media_type FROM _bioextract.source_file"
        ).fetchone() == ("application/vnd.apache.parquet",)

    hive = tmp_path / "hive" / "tax_id=9606"
    hive.mkdir(parents=True)
    frame.drop("tax_id").write_parquet(hive / "part.parquet")
    (tmp_path / "hive" / "README.txt").write_text("ignored", encoding="utf-8")
    assert (
        UniProtDatabase.from_idmapping(tmp_path / "hive")
        .read_mapping(taxon_ids=["9606"])
        .height
        == 1
    )


def test_idmapping_taxon_validation_and_all_taxa_opt_in(tmp_path: Path) -> None:
    database = UniProtDatabase.from_idmapping(
        _write_idmapping(tmp_path / "mapping.tab.gz")
    )
    with pytest.raises(ValueError, match="tax_id"):
        database.scan_mapping(taxon_ids=["9606", ""])
    with pytest.raises(ValueError, match="allow_all_taxa"):
        database.write_duckdb(tmp_path / "blocked.duckdb")
    assert database.write_duckdb(
        tmp_path / "all.duckdb", allow_all_taxa=True
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
    assert matches.extract_proteins()["primary_accession"].to_list() == ["P12345"]
    assert matches.extract_proteins()["protein_existence"].to_list() == [
        "1: Evidence at protein level"
    ]
    assert "group_id" not in matches.extract_proteins().schema
    assert matches.extract_unmatched_ids()["input_id"].to_list() == [
        "TEST",
        "missing",
    ]
    assert database.select_ids(["TEST"], namespace="gene_name").extract_gene_names()[
        "name"
    ].to_list() == ["TEST", "ALT", "LOC1", "ORF1", "TEST2", "ALT2"]
    assert database.select_ids(["P12345"], namespace="uniprot").extract_protein_names()[
        "name"
    ].to_list() == [
        "Test protein",
        "First alternative",
        "Second alternative",
    ]
    assert database.select_ids(["1234"], namespace="gene_id").extract_cross_references(
        databases=["GeneID"]
    )["external_id"].to_list() == ["1234"]
    refseq = database.select_ids(
        ["NP_000001.1"], namespace="refseq"
    ).extract_cross_references(databases=["RefSeq"])
    assert refseq.select("external_id", "isoform_id").to_dicts() == [
        {"external_id": "NP_000001.1", "isoform_id": "P12345-2"}
    ]
    con_server = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_cross_references(databases=["ConoServer"])
    assert con_server.select("external_id", "properties", "isoform_id").to_dicts() == [
        {
            "external_id": "1570",
            "properties": "GIIIA [R13A]",
            "isoform_id": None,
        }
    ]
    sequences = database.select_ids(["P12345"], namespace="uniprot").extract_sequences(
        sequence_type="all"
    )
    assert sequences.schema["crc64"] == pl.String
    assert sequences.select("sequence_type", "crc64").to_dicts() == [
        {"sequence_type": "canonical", "crc64": "81FC3551E879CB1A"},
        {"sequence_type": "isoform", "crc64": None},
    ]
    selection = database.select_ids(["P12345"], namespace="uniprot")
    assert selection.extract_ec_numbers()["ec_number"].to_list() == ["1.2.3.4"]
    assert selection.extract_go_annotations()["go_id"].to_list() == ["GO:0003677"]
    assert selection.extract_comments(comment_types=["FUNCTION"])[
        ["comment_type", "comment_text"]
    ].to_dicts() == [
        {
            "comment_type": "FUNCTION",
            "comment_text": "Demonstrates a mem-brane fixture.",
        }
    ]
    assert selection.extract_subcellular_locations()["location"].to_list() == [
        "Nucleus"
    ]
    assert selection.extract_keywords()["keyword"].to_list() == [
        "Reference proteome",
        "Nucleus",
    ]
    grouped = database.select_groups(
        {"case": ["P12345"], "control": ["missing"]},
        namespace="uniprot",
    )
    assert grouped.extract_proteins()["group_id"].to_list() == ["case"]
    assert grouped.extract_unmatched_ids().select(
        "group_id", "input_id"
    ).to_dicts() == [{"group_id": "control", "input_id": "missing"}]
    isoforms_frame = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_isoforms()
    assert isoforms_frame.select(
        "isoform_id", "isoform_order", "sequence_status", "sequence_id"
    ).to_dicts() == [
        {
            "isoform_id": "P12345-1",
            "isoform_order": 1,
            "sequence_status": "Displayed",
            "sequence_id": "P12345:canonical",
        },
        {
            "isoform_id": "P12345-2",
            "isoform_order": 2,
            "sequence_status": "Alternative",
            "sequence_id": "P12345-2",
        },
    ]
    assert selection.extract_isoform_identifiers().select(
        "isoform_id", "identifier", "identifier_order", "is_main"
    ).to_dicts() == [
        {
            "isoform_id": "P12345-1",
            "identifier": "P12345-1",
            "identifier_order": 1,
            "is_main": True,
        },
        {
            "isoform_id": "P12345-2",
            "identifier": "P12345-2",
            "identifier_order": 1,
            "is_main": True,
        },
    ]
    assert database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_isoform_variations().select(
        "isoform_id", "variation_id", "variation_order"
    ).to_dicts() == [
        {
            "isoform_id": "P12345-2",
            "variation_id": "VSP_000001",
            "variation_order": 1,
        }
    ]
    variation = database.select_ids(
        ["P12345"], namespace="uniprot"
    ).extract_sequence_variations()
    assert variation.select("variation_id", "note").to_dicts() == [
        {"variation_id": "VSP_000001", "note": "AA -> GG (in isoform 2)"}
    ]
    assert database.select_ids(
        ["P12345"], namespace="uniprot", taxon_ids=["10090"]
    ).extract_proteins().columns == [
        "input_id",
        "input_namespace",
        "primary_accession",
        "entry_name",
        "is_reviewed",
        "taxon_id",
        "protein_existence",
        "sequence_length",
        "molecular_weight",
        "sequence_version",
        "entry_version",
    ]
    assert database.select_ids(
        [], namespace="uniprot"
    ).extract_accessions().columns == [
        "input_id",
        "input_namespace",
        "primary_accession",
        "accession",
        "accession_order",
        "is_primary",
    ]
    with database.connect() as connection:
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
        assert metadata["bioextract.metadata_schema_version"] == "1"
        assert (
            metadata["bioextract.resource_schema_version"]
            == "uniprot-knowledgebase-duckdb-v1"
        )
        assert metadata["bioextract.source_schema_profile"] == "uniprotkb-flat-file-v1"
        assert (
            metadata["bioextract.molecular_weight_validation_model"]
            == "compatible-current-and-legacy"
        )
        with pytest.raises(duckdb.Error):
            connection.execute("CREATE TABLE forbidden(value INTEGER)")


@pytest.mark.parametrize(
    ("role", "records"),
    [
        ("canonical_sequences", {"P12345": "ACDEFGHIKL"}),
        ("isoform_sequences", {"P12345-2": "ACGGFGHIKL"}),
    ],
)
def test_fasta_plain_and_gzip_inputs_are_equivalent(
    tmp_path: Path, role: str, records: dict[str, str]
) -> None:
    entries = _write_dat(tmp_path / "entries.dat.gz")
    plain = _write_fasta(tmp_path / f"{role}.fasta", records, compressed=False)
    compressed = _write_fasta(tmp_path / f"{role}.compressed", records)
    plain_publication = tmp_path / f"{role}-plain.duckdb"
    compressed_publication = tmp_path / f"{role}-compressed.duckdb"

    UniProtDatabase.from_knowledgebase(entries=entries, **{role: plain}).write_duckdb(
        plain_publication
    )
    UniProtDatabase.from_knowledgebase(
        entries=entries, **{role: compressed}
    ).write_duckdb(compressed_publication)

    for table in ("protein_sequence", "protein_isoform"):
        with (
            duckdb.connect(str(plain_publication), read_only=True) as plain_connection,
            duckdb.connect(
                str(compressed_publication), read_only=True
            ) as compressed_connection,
        ):
            plain_rows = plain_connection.execute(
                f"SELECT * FROM {table} ORDER BY ALL"
            ).fetchall()
            compressed_rows = compressed_connection.execute(
                f"SELECT * FROM {table} ORDER BY ALL"
            ).fetchall()
        assert plain_rows == compressed_rows


@pytest.mark.parametrize(
    ("role", "role_label"),
    [("canonical_sequences", "canonical"), ("isoform_sequences", "isoform")],
)
def test_corrupt_gzip_fasta_error_includes_role_and_path(
    tmp_path: Path, role: str, role_label: str
) -> None:
    entries = _write_dat(tmp_path / "entries.dat.gz")
    corrupt = tmp_path / f"{role}.fasta.gz"
    corrupt.write_bytes(b"\x1f\x8b\x08\x00not-a-valid-gzip")

    with pytest.raises(ValueError) as raised:
        UniProtDatabase.from_knowledgebase(
            entries=entries, **{role: corrupt}
        ).write_duckdb(tmp_path / f"{role}-corrupt.duckdb")

    message = str(raised.value)
    assert f"UniProt {role_label} FASTA input" in message
    assert "invalid gzip stream" in message
    assert str(corrupt) in message


def test_shared_secondary_accession_selects_every_canonical_match(
    tmp_path: Path,
) -> None:
    entries = tmp_path / "shared-secondary-accession.dat"
    entries.write_text(
        """ID   FIRST_HUMAN Reviewed; 3 AA.
AC   P11111; Q99999;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
     ACD
//
ID   SECOND_MOUSE Reviewed; 3 AA.
AC   P22222; Q99999;
OX   NCBI_TaxID=10090;
SQ   SEQUENCE   3 AA;  365 MW;  69CB1DB000000000 CRC64;
     AEF
//
""",
        encoding="utf-8",
    )
    publication = tmp_path / "shared-secondary.duckdb"
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(publication)
    database = UniProtDatabase.from_duckdb(publication)

    assert database.select_ids(
        ["Q99999"], namespace="uniprot"
    ).extract_proteins().select(
        "input_id", "primary_accession", "taxon_id"
    ).to_dicts() == [
        {"input_id": "Q99999", "primary_accession": "P11111", "taxon_id": "9606"},
        {"input_id": "Q99999", "primary_accession": "P22222", "taxon_id": "10090"},
    ]
    assert database.select_ids(
        ["Q99999"], namespace="uniprot", taxon_ids=["10090"]
    ).extract_proteins()["primary_accession"].to_list() == ["P22222"]
    assert database.select_groups(
        {"demerged": ["Q99999"]}, namespace="uniprot"
    ).extract_proteins().select("group_id", "primary_accession").to_dicts() == [
        {"group_id": "demerged", "primary_accession": "P11111"},
        {"group_id": "demerged", "primary_accession": "P22222"},
    ]


def test_group_selection_resolves_each_identifier_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    publication = tmp_path / "uniprot.duckdb"
    UniProtDatabase.from_knowledgebase(
        entries=_write_dat(tmp_path / "entries.dat.gz")
    ).write_duckdb(publication)
    database = UniProtDatabase.from_duckdb(publication)
    query_calls: list[tuple[str, ...]] = []
    original_query = (
        uniprot_query.UniProtSelection._query_identifier_matches  # pyright: ignore[reportPrivateUsage]
    )

    def counted_query(
        selection: uniprot_query.UniProtSelection,
    ) -> pl.DataFrame:
        query_calls.append(selection.input_ids)
        return original_query(selection)

    monkeypatch.setattr(
        uniprot_query.UniProtSelection,
        "_query_identifier_matches",
        counted_query,
    )
    selection = database.select_groups(
        {
            "case": [" P12345 ", "P12345", "missing"],
            "control": ["P12345", "missing"],
            "empty": [],
        },
        namespace="uniprot",
    )

    assert selection.input_ids == ("P12345", "missing")
    assert selection.group_ids == ("case", "control", "empty")
    assert selection.extract_proteins().select("group_id", "input_id").to_dicts() == [
        {"group_id": "case", "input_id": "P12345"},
        {"group_id": "control", "input_id": "P12345"},
    ]
    selection.extract_accessions()
    assert selection.extract_unmatched_ids().select(
        "group_id", "input_id"
    ).to_dicts() == [
        {"group_id": "case", "input_id": "missing"},
        {"group_id": "control", "input_id": "missing"},
    ]
    assert query_calls == [("P12345", "missing")]


def test_selection_and_fasta_inputs_reject_empty_values(tmp_path: Path) -> None:
    path = tmp_path / "uniprot.duckdb"
    entries = _write_dat(tmp_path / "entries.dat.gz")
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(path)
    database = UniProtDatabase.from_duckdb(path)
    with pytest.raises(ValueError, match="tax_id"):
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

    leading_whitespace = tmp_path / "leading-whitespace.fasta"
    leading_whitespace.write_text(
        ">sp|P12345|TEST_HUMAN\n ACDEFGHIKL\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="invalid sequence characters"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, canonical_sequences=leading_whitespace
        ).write_duckdb(tmp_path / "leading-whitespace.duckdb")

    trailing_whitespace = tmp_path / "trailing-whitespace.fasta"
    trailing_whitespace.write_text(">sp|P12345-2|TEST_HUMAN\nACD \n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid sequence characters"):
        UniProtDatabase.from_knowledgebase(
            entries=entries, isoform_sequences=trailing_whitespace
        ).write_duckdb(tmp_path / "trailing-whitespace.duckdb")

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


def test_knowledgebase_accepts_real_q6t412_rounding_vector(
    tmp_path: Path,
) -> None:
    sequence = (
        "MTEVQPPPAQSTVATADTPSLAPDTTLETSTSTELAPITTEQTIITTNAEGKKVKKIIRR"
        "KRRPARPQVDPATFKTDT PAPT GTSFNIWYNKWSGGDREDKYLSQTAAQGRCNVARDSGY"
        "TKADKTPGSYFCLFFARGICPKGVDCEYLHRLPTVTDIFPSNIDCFGRDKHSDYRDDMGG"
        "VGSFQRQNRTLYIGRIHVTDDIEEIVARHFQEWGQIERTRVLTARGVAFVTYMNEANSQF"
        "AKEAMAHQSLDHNEILNVRWATVDPNPQAAKREAHRIEEQAAEAIRKALPAAYVAELEGR"
        "DPEAKKRRKIEGSFGLQGYEAPDDVWYAKEKAEWEAAKEIEAAGGAAXPRQMIESGEDAH"
        "AHEADCAAMQVAPSGQHSQGNGIFSTSTLAALRGYTAAPAKPKVAPVAGPLVGYGSDDDSD"
    ).replace(" ", "")
    assert len(sequence) == 421

    entries = tmp_path / "q6t412.dat"
    entries.write_text(
        "ID   Q6T412_FIXTURE Reviewed; 421 AA.\n"
        "AC   Q6T412;\n"
        "OX   NCBI_TaxID=9606;\n"
        "SQ   SEQUENCE   421 AA;  46189 MW;  EE7D1FA88E010B94 CRC64;\n"
        f"     {sequence}\n"
        "//\n",
        encoding="utf-8",
    )
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
        tmp_path / "q6t412.duckdb"
    )


def test_dat_molecular_weight_must_match_sequence(tmp_path: Path) -> None:
    entries = tmp_path / "bad-molecular-weight.dat"
    entries.write_text(
        """ID   TEST_HUMAN Reviewed; 3 AA.
AC   P12345;
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   3 AA;  308 MW;  6AAEBDB000000000 CRC64;
     ACD
//
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="molecular weight mismatch"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "bad-molecular-weight.duckdb"
        )


@pytest.mark.parametrize("molecular_weight", [5719, 5720])
def test_knowledgebase_accepts_one_consistent_sec_weight_model(
    tmp_path: Path, molecular_weight: int
) -> None:
    entries = tmp_path / f"sec-{molecular_weight}.dat"
    entries.write_text(
        "ID   SEC_FIXTURE Reviewed; 38 AA.\n"
        f"AC   P{molecular_weight};\n"
        "OX   NCBI_TaxID=9606;\n"
        f"SQ   SEQUENCE   38 AA;  {molecular_weight} MW;  "
        "3100000000007707 CRC64;\n"
        f"     {'U' * 38}\n"
        "//\n",
        encoding="utf-8",
    )
    UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
        tmp_path / f"sec-{molecular_weight}.duckdb"
    )


def test_knowledgebase_rejects_conflicting_sec_weight_models(tmp_path: Path) -> None:
    sequence = "U" * 38
    entries = tmp_path / "conflicting-sec.dat"
    entries.write_text(
        "ID   CURRENT_SEC Reviewed; 38 AA.\n"
        "AC   P5720;\n"
        "OX   NCBI_TaxID=9606;\n"
        "SQ   SEQUENCE   38 AA;  5720 MW;  3100000000007707 CRC64;\n"
        f"     {sequence}\n"
        "//\n"
        "ID   LEGACY_SEC Reviewed; 38 AA.\n"
        "AC   P5719;\n"
        "OX   NCBI_TaxID=9606;\n"
        "SQ   SEQUENCE   38 AA;  5719 MW;  3100000000007707 CRC64;\n"
        f"     {sequence}\n"
        "//\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="conflicting molecular-weight models"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "conflicting-sec.duckdb"
        )


@pytest.mark.parametrize(
    ("alternative_products", "feature", "message"),
    [
        (
            "CC       Name=Bad; IsoId=P12345-2; Sequence=Bogus;\n",
            "",
            "alternative-products Sequence",
        ),
        (
            "CC       Name=Bad; IsoId=P12345-2; Sequence=VSP_000001;\n",
            'FT   VAR_SEQ         2\nFT                   /note="C -> G"\n',
            "VAR_SEQ feature is missing /id",
        ),
        (
            "CC       Name=Bad; Sequence=Displayed;\n",
            "",
            "Name block is missing IsoId",
        ),
        (
            "CC       Name=Bad; IsoId=P12345-2;\n",
            "",
            "Name block is missing Sequence",
        ),
        (
            "CC       Name=Bad; IsoId=P12345-2; Sequence=VSP_000001;\n",
            'FT   VAR_SEQ         2\nFT                   /id="VSP_1"\n',
            "Invalid UniProt VAR_SEQ /id",
        ),
        (
            "CC       Name=Bad; IsoId=P12345-2; Sequence=VSP_000001;\n",
            (
                'FT   VAR_SEQ         2\nFT                   /id="VSP_000001"\n'
                'FT                   /id="VSP_000001"\n'
            ),
            "Duplicate UniProt VAR_SEQ /id",
        ),
    ],
)
def test_knowledgebase_rejects_invalid_isoform_sequence_semantics(
    tmp_path: Path,
    alternative_products: str,
    feature: str,
    message: str,
) -> None:
    entries = tmp_path / "invalid-isoform.dat"
    entries.write_text(
        (
            "ID   TEST_HUMAN Reviewed; 3 AA.\n"
            "AC   P12345;\n"
            "OX   NCBI_TaxID=9606;\n"
            "CC   -!- ALTERNATIVE PRODUCTS:\n"
            f"{alternative_products}"
            f"{feature}"
            "SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;\n"
            "     ACD\n"
            "//\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=message):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "invalid-isoform.duckdb"
        )


@pytest.mark.parametrize(
    "sq_line",
    [
        "SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDB0 CRC64;",
        "SQ   SEQUENCE   3 AA;  300 MW;  6aaebdb000000000 CRC64;",
        "SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDBG000000000 CRC64;",
        "SQ   SEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64; trailing",
        "SQ\tSEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64;",
        "SQ   junk SEQUENCE   3 AA;  300 MW;  6AAEBDB000000000 CRC64;",
    ],
)
def test_dat_sq_line_requires_exact_crc64_grammar(tmp_path: Path, sq_line: str) -> None:
    entries = tmp_path / "invalid-sq.dat"
    entries.write_text(
        "\n".join(
            [
                "ID   TEST_HUMAN Reviewed; 3 AA.",
                "AC   P12345;",
                "OX   NCBI_TaxID=9606;",
                sq_line,
                "     ACD",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Invalid UniProt SQ line"):
        UniProtDatabase.from_knowledgebase(entries=entries).write_duckdb(
            tmp_path / "invalid-sq.duckdb"
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
SQ   SEQUENCE   3 AA;  307 MW;  6AAEBDB000000000 CRC64;
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
SQ   SEQUENCE   3 AA;  365 MW;  69CB1DB000000000 CRC64;
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
        "primary_accession"
    ].to_list() == ["P16375", "P16376"]
    assert database.select_ids(
        ["P16375"], namespace="uniprot"
    ).extract_isoform_variations()["variation_id"].to_list() == ["VSP_013348"]
    assert (
        database.select_ids(["P16376"], namespace="uniprot")
        .extract_isoform_variations()
        .is_empty()
    )
