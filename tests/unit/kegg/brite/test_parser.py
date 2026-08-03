from bioextract.kegg.brite.parse import parse_entry_and_ko, parse_pathway_level3


def test_parse_pathway_level3_supports_path_br_and_plain_forms() -> None:
    level3_path = parse_pathway_level3(
        "00010 Glycolysis / Gluconeogenesis [PATH:tcar00010]"
    )
    assert level3_path.id == "00010"
    assert level3_path.kegg_id == "tcar00010"
    assert level3_path.name == "Glycolysis / Gluconeogenesis"

    level3_plain = parse_pathway_level3("00566 Sulfoquinovose metabolism")
    assert level3_plain.id == "00566"
    assert level3_plain.kegg_id is None
    assert level3_plain.name == "Sulfoquinovose metabolism"

    level3_br = parse_pathway_level3("01001 Protein kinases [BR:tcar01001]")
    assert level3_br.id == "01001"
    assert level3_br.kegg_id == "tcar01001"
    assert level3_br.name == "Protein kinases"


def test_parse_entry_and_ko_supports_entry_only_leaf() -> None:
    leaf = parse_entry_and_ko("U0034_01605 oxidoreductase")
    assert leaf.entry.id == "U0034_01605"
    assert leaf.entry.name == "oxidoreductase"
    assert leaf.ko is None
