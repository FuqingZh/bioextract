from bioextract.kegg.metabolic.core import (
    ModuleDefinitionParser,
    ast_rows,
    parse_equation,
)


def test_parse_equation_ignores_participant_suffix_qualifiers() -> None:
    assert parse_equation(
        "R00379",
        "C00039(n+1) + G10477(n) <=> 2 C01330(side 1) + C00001",
    ) == [
        {
            "reaction_id": "R00379",
            "side": "left",
            "position": 1,
            "participant_namespace": "kegg_compound",
            "participant_id": "C00039",
            "coefficient_text": "1",
            "coefficient_numeric": 1.0,
        },
        {
            "reaction_id": "R00379",
            "side": "left",
            "position": 2,
            "participant_namespace": "kegg_glycan",
            "participant_id": "G10477",
            "coefficient_text": "1",
            "coefficient_numeric": 1.0,
        },
        {
            "reaction_id": "R00379",
            "side": "right",
            "position": 1,
            "participant_namespace": "kegg_compound",
            "participant_id": "C01330",
            "coefficient_text": "2",
            "coefficient_numeric": 2.0,
        },
        {
            "reaction_id": "R00379",
            "side": "right",
            "position": 2,
            "participant_namespace": "kegg_compound",
            "participant_id": "C00001",
            "coefficient_text": "1",
            "coefficient_numeric": 1.0,
        },
    ]


def test_module_parser_treats_double_dash_as_optional_placeholder() -> None:
    rows = ast_rows(
        "M00076",
        ModuleDefinitionParser("-- K01136 K01217").parse(),
    )

    assert [
        (row["node_kind"], row["member_namespace"], row["member_id"]) for row in rows
    ] == [
        ("sequence", None, None),
        ("optional", None, None),
        ("identifier", "ko", "K01136"),
        ("identifier", "ko", "K01217"),
    ]
