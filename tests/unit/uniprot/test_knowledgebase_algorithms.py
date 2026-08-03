from __future__ import annotations

import sqlite3

import bioextract.uniprot._knowledgebase as knowledgebase


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


def test_uniprot_molecular_weight_matches_real_q6gzx4_vector() -> None:
    sequence = (
        "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS"
        "EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD"
        "AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL"
        "EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD"
        "SFRKIYTDLGWKFTPL"
    )
    assert len(sequence) == 256
    assert (
        knowledgebase._calculate_molecular_weight(  # pyright: ignore[reportPrivateUsage]
            sequence
        )
        == 29735
    )


def test_uniprot_molecular_weight_matches_real_q6t412_rounding_vector() -> None:
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
    assert (
        knowledgebase._calculate_molecular_weight(  # pyright: ignore[reportPrivateUsage]
            sequence
        )
        == 46189
    )
