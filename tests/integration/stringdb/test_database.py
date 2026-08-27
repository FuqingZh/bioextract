from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

import bioextract.stringdb.stringdb as stringdb_module
from bioextract.errors import CapabilityError, IntegrityError
from bioextract.stringdb import STRINGDatabase


def _write_text_or_gzip(path: Path, content: str, *, should_gzip: bool) -> None:
    if should_gzip:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return

    path.write_text(content, encoding="utf-8")


def _write_demo_string_files(
    *,
    aliases: Path,
    links: Path,
    should_gzip: bool = False,
) -> None:
    _write_text_or_gzip(
        aliases,
        "\n".join(
            [
                "string_protein_id\talias\tsource",
                "9606.ENSP0001\tP04637\tUniProt_AC",
                "9606.ENSP0001\tTP53\tUniProt_GN_Name",
                "9606.ENSP0002\tEGFR\tUniProt_GN_Name",
                "9606.ENSP0003\tCDK2\tUniProt_GN_Name",
                "9606.ENSP9999\tP04637\tUniProt_GN_Synonyms",
            ]
        )
        + "\n",
        should_gzip=should_gzip,
    )
    _write_text_or_gzip(
        links,
        "\n".join(
            [
                "protein1 protein2 combined_score",
                "9606.ENSP0001 9606.ENSP0002 700",
                "9606.ENSP0002 9606.ENSP0001 700",
                "9606.ENSP0001 9606.ENSP0001 999",
                "9606.ENSP0001 9606.ENSP0003 450",
                "9606.ENSP0002 9606.ENSP9999 800",
            ]
        )
        + "\n",
        should_gzip=should_gzip,
    )


def test_extract_string_mapping_accepts_hash_string_id_header(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.tsv"
    file_links = tmp_path / "links.tsv"
    _write_text_or_gzip(
        file_aliases,
        "\n".join(
            [
                "#string_protein_id\talias\tsource",
                "9606.ENSP0001\tP04637\tUniProt_AC",
                "9606.ENSP0002\tEGFR\tUniProt_GN_Synonyms",
                "9606.ENSP0002\tEGFR\tUniProt_GN_Name",
                "9606.ENSP0003\tCDK2\tKEGG_GENEID",
            ]
        )
        + "\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 700\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(
        aliases=file_aliases,
        links=file_links,
    ).select_ids([" sp|P04637|P53_HUMAN ", " EGFR ", "  "])

    df_result = selection.mappings().collect()

    assert df_result.columns == [
        "input_id",
        "#string_protein_id",
        "alias",
        "source",
    ]
    assert df_result.to_dicts() == [
        {
            "input_id": "EGFR",
            "#string_protein_id": "9606.ENSP0002",
            "alias": "EGFR",
            "source": "UniProt_GN_Name",
        },
        {
            "input_id": "P04637",
            "#string_protein_id": "9606.ENSP0001",
            "alias": "P04637",
            "source": "UniProt_AC",
        },
    ]


def test_stringdb_single_query_minimal_round_trip(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_min_combined_score(400)
    )

    assert selection.edges().collect().columns == [
        "string_id_a",
        "string_id_b",
        "combined_score",
    ]
    assert selection.mappings().collect().columns == [
        "input_id",
        "string_protein_id",
        "alias",
        "source",
    ]
    assert selection.unmatched_ids().collect().columns == ["input_id"]
    assert selection.edges().collect().to_dicts() == [
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
    ]


def test_selection_exposes_group_mode(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)
    db = STRINGDatabase.from_files(aliases=file_aliases, links=file_links)

    assert db.select_ids(["TP53"]).is_grouped is False
    assert db.select_groups({"G1": ["TP53"]}).is_grouped is True


def test_string_alias_pipe_grammar_does_not_rewrite_direct_string_ids(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)
    database = STRINGDatabase.from_files(
        aliases=file_aliases,
        links=file_links,
    )

    for invalid in ("db|P04637|P53", "sp||P53", "sp|P04637|"):
        with pytest.raises(ValueError, match="UniProt pipe-form"):
            database.select_ids([invalid], namespace="alias")
        with pytest.raises(ValueError, match="UniProt pipe-form"):
            database.select_groups({"case": [invalid]}, namespace="alias")

    exact_string = database.select_ids(
        [" sp|P04637|P53_HUMAN "],
        namespace="string",
    )
    assert exact_string.unmatched_ids().collect().to_dicts() == [
        {"input_id": "sp|P04637|P53_HUMAN"}
    ]


def test_stringdb_group_query_minimal_round_trip(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_groups({"G1": ["TP53", "EGFR"], "G2": ["P04637", "CDK2"]})
        .with_min_combined_score(400)
    )

    df_edges = selection.edges().collect()
    df_mapping = selection.mappings().collect()
    df_unmapped = selection.unmatched_ids().collect()

    assert df_edges.columns == [
        "group_id",
        "string_id_a",
        "string_id_b",
        "combined_score",
    ]
    assert df_mapping.columns == [
        "group_id",
        "input_id",
        "string_protein_id",
        "alias",
        "source",
    ]
    assert df_unmapped.columns == ["group_id", "input_id"]
    assert df_edges.to_dicts() == [
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
        {
            "group_id": "G2",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "group_id": "G2",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0003",
            "combined_score": 450,
        },
    ]


def test_extract_string_mapping_accepts_plain_and_gzip_inputs(tmp_path: Path) -> None:
    for should_gzip in (False, True):
        suffix = ".txt.gz" if should_gzip else ".txt"
        file_aliases = tmp_path / f"aliases{suffix}"
        file_links = tmp_path / f"links{suffix}"
        _write_text_or_gzip(
            file_aliases,
            "\n".join(
                [
                    "string_protein_id\talias\tsource",
                    "9606.ENSP0001\tTP53\tUniProt_GN_Name",
                    "9606.ENSP0002\tEGFR\tUniProt_GN_Name",
                ]
            )
            + "\n",
            should_gzip=should_gzip,
        )
        _write_text_or_gzip(
            file_links,
            "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 500\n",
            should_gzip=should_gzip,
        )

        selection = STRINGDatabase.from_files(
            aliases=file_aliases,
            links=file_links,
        ).select_ids(["TP53", "EGFR"])

        df_result = selection.mappings().collect()
        assert df_result.to_dicts() == [
            {
                "input_id": "EGFR",
                "string_protein_id": "9606.ENSP0002",
                "alias": "EGFR",
                "source": "UniProt_GN_Name",
            },
            {
                "input_id": "TP53",
                "string_protein_id": "9606.ENSP0001",
                "alias": "TP53",
                "source": "UniProt_GN_Name",
            },
        ]


def test_extract_core_outputs_return_expected_frames(tmp_path: Path) -> None:
    for should_gzip in (False, True):
        suffix = ".txt.gz" if should_gzip else ".txt"
        file_aliases = tmp_path / f"aliases{suffix}"
        file_links = tmp_path / f"links{suffix}"
        _write_demo_string_files(
            aliases=file_aliases,
            links=file_links,
            should_gzip=should_gzip,
        )

        selection = (
            STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
            .select_ids(["sp|P04637|P53_HUMAN", "EGFR", "CDK2", "MISSING"])
            .with_min_combined_score(300)
        )
        df_edges = selection.edges().collect()
        df_mapping = selection.mappings().collect()
        df_unmapped = selection.unmatched_ids().collect()

        assert df_edges.columns == ["string_id_a", "string_id_b", "combined_score"]
        assert df_mapping.columns == [
            "input_id",
            "string_protein_id",
            "alias",
            "source",
        ]
        assert df_unmapped.columns == ["input_id"]

        assert df_edges.to_dicts() == [
            {
                "string_id_a": "9606.ENSP0001",
                "string_id_b": "9606.ENSP0001",
                "combined_score": 999,
            },
            {
                "string_id_a": "9606.ENSP0001",
                "string_id_b": "9606.ENSP0002",
                "combined_score": 700,
            },
            {
                "string_id_a": "9606.ENSP0001",
                "string_id_b": "9606.ENSP0003",
                "combined_score": 450,
            },
            {
                "string_id_a": "9606.ENSP0002",
                "string_id_b": "9606.ENSP9999",
                "combined_score": 800,
            },
        ]
        assert df_mapping.to_dicts() == [
            {
                "input_id": "CDK2",
                "string_protein_id": "9606.ENSP0003",
                "alias": "CDK2",
                "source": "UniProt_GN_Name",
            },
            {
                "input_id": "EGFR",
                "string_protein_id": "9606.ENSP0002",
                "alias": "EGFR",
                "source": "UniProt_GN_Name",
            },
            {
                "input_id": "P04637",
                "string_protein_id": "9606.ENSP0001",
                "alias": "P04637",
                "source": "UniProt_AC",
            },
            {
                "input_id": "P04637",
                "string_protein_id": "9606.ENSP9999",
                "alias": "P04637",
                "source": "UniProt_GN_Synonyms",
            },
        ]
        assert df_unmapped.to_dicts() == [{"input_id": "MISSING"}]


def test_extract_core_outputs_handle_empty_input_set(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection = STRINGDatabase.from_files(
        aliases=file_aliases, links=file_links
    ).select_ids([])
    df_edges = selection.edges().collect()
    df_mapping = selection.mappings().collect()
    df_unmapped = selection.unmatched_ids().collect()

    assert df_edges.schema == {
        "string_id_a": pl.String,
        "string_id_b": pl.String,
        "combined_score": pl.Int64,
    }
    assert df_mapping.schema == {
        "input_id": pl.String,
        "string_protein_id": pl.String,
        "alias": pl.String,
        "source": pl.String,
    }
    assert df_unmapped.schema == {"input_id": pl.String}


def test_extract_string_mapping_honors_custom_source_rank_map(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_aliases,
        "\n".join(
            [
                "string_protein_id\talias\tsource",
                "9606.ENSP0001\tTP53\tUniProt_GN_Name",
                "9606.ENSP0001\tTP53\tCustomSource",
            ]
        )
        + "\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 500\n",
        should_gzip=False,
    )

    df_default = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_ids(["TP53"])
        .mappings()
        .collect()
    )
    df_custom = (
        STRINGDatabase.from_files(
            aliases=file_aliases,
            links=file_links,
            rank_by_source={
                "CustomSource": 1,
                **STRINGDatabase.DEFAULT_SOURCE_RANK_MAP,
            },
        )
        .select_ids(["TP53"])
        .mappings()
        .collect()
    )

    assert df_default.to_dicts() == [
        {
            "input_id": "TP53",
            "string_protein_id": "9606.ENSP0001",
            "alias": "TP53",
            "source": "UniProt_GN_Name",
        }
    ]
    assert df_custom.to_dicts() == [
        {
            "input_id": "TP53",
            "string_protein_id": "9606.ENSP0001",
            "alias": "TP53",
            "source": "CustomSource",
        }
    ]


def test_selection_reuses_cached_frames(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_min_combined_score(400)
    )

    df_protein_map_first = selection.mappings().collect()
    df_protein_map_second = selection.mappings().collect()
    assert df_protein_map_first.equals(df_protein_map_second)

    df_edges_first = selection.edges().collect()
    df_edges_second = selection.edges().collect()
    assert df_edges_first.equals(df_edges_second)


def test_score_filter_reuses_mapping_result_but_recomputes_edges(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection_low = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_min_combined_score(300)
    )
    selection_low_map = selection_low.mappings().collect()
    selection_low_unmapped = selection_low.unmatched_ids().collect()
    selection_low_edges = selection_low.edges().collect()

    selection_high = selection_low.with_min_combined_score(600)

    assert selection_high.mappings().collect().equals(selection_low_map)
    assert selection_high.unmatched_ids().collect().equals(selection_low_unmapped)
    # This fixture's induced network contains only scores >= 600, so the
    # threshold change preserves the edge result while still producing a new
    # lazy plan.
    assert selection_high.edges().collect().equals(selection_low_edges)
    assert selection_low.edges().collect().to_dicts() == [
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
    ]
    assert selection_high.edges().collect().to_dicts() == [
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
    ]


def test_from_files_rejects_missing_files(tmp_path: Path) -> None:
    file_aliases = tmp_path / "missing_aliases.txt"
    file_links = tmp_path / "missing_links.txt"

    with pytest.raises(FileNotFoundError, match="aliases file not found"):
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)

    file_aliases.write_text("", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="links file not found"):
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)


def test_from_files_rejects_empty_release_version(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    with pytest.raises(ValueError, match="release_version"):
        STRINGDatabase.from_files(
            aliases=file_aliases,
            links=file_links,
            release_version="",
        )


def test_select_groups_rejects_invalid_group_shape(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    db = STRINGDatabase.from_files(aliases=file_aliases, links=file_links)

    with pytest.raises(ValueError, match="non-empty"):
        db.select_groups({"  ": ["TP53"]})

    with pytest.raises(ValueError, match="unique after normalization"):
        db.select_groups({"A": ["TP53"], " A ": ["EGFR"]})


def test_extract_string_mapping_rejects_aliases_missing_required_columns(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_aliases,
        "bad_column\talias\tsource\nfoo\tTP53\tUniProt_GN_Name\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 500\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(
        aliases=file_aliases,
        links=file_links,
    ).select_ids(["TP53"])

    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="aliases file"):
        selection.mappings().collect()


def test_extract_edges_rejects_links_missing_required_columns(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_aliases,
        "string_protein_id\talias\tsource\n9606.ENSP0001\tTP53\tUniProt_GN_Name\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2\n9606.ENSP0001 9606.ENSP0002\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(
        aliases=file_aliases,
        links=file_links,
    ).select_ids(["TP53"])

    with pytest.raises((ValueError, pl.exceptions.ComputeError), match="links file"):
        selection.edges().collect()


def test_extract_string_mapping_works_without_links_file(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    _write_text_or_gzip(
        file_aliases,
        "string_protein_id\talias\tsource\n9606.ENSP0001\tTP53\tUniProt_GN_Name\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(aliases=file_aliases).select_groups(
        {"G1": ["TP53"], "G2": ["MISSING"]}
    )

    assert selection.mappings().collect().to_dicts() == [
        {
            "group_id": "G1",
            "input_id": "TP53",
            "string_protein_id": "9606.ENSP0001",
            "alias": "TP53",
            "source": "UniProt_GN_Name",
        }
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "G2", "input_id": "MISSING"}
    ]


def test_extract_edges_rejects_missing_links_file(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    _write_text_or_gzip(
        file_aliases,
        "string_protein_id\talias\tsource\n9606.ENSP0001\tTP53\tUniProt_GN_Name\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(aliases=file_aliases).select_ids(["TP53"])
    with pytest.raises(
        CapabilityError, match="links source is absent from this snapshot"
    ):
        selection.dataset.scan_links()

    with pytest.raises(
        (CapabilityError, pl.exceptions.ComputeError),
        match="links source is absent from this snapshot",
    ):
        selection.edges().collect()


def test_extract_string_mapping_rejects_missing_aliases_file(tmp_path: Path) -> None:
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 500\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(links=file_links).select_ids(["TP53"])
    with pytest.raises(
        CapabilityError, match="aliases source is absent from this snapshot"
    ):
        selection.dataset.scan_aliases()

    with pytest.raises(
        (CapabilityError, pl.exceptions.ComputeError),
        match="aliases source is absent from this snapshot",
    ):
        selection.mappings().collect()


def test_group_selection_extracts_flat_outputs(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    group_selection = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_groups(
            {
                "G1": ["sp|P04637|P53_HUMAN", "EGFR", "MISSING"],
                "G2": ["P04637", "CDK2"],
                "G3": [],
            }
        )
        .with_min_combined_score(300)
    )

    df_group_mapping = group_selection.mappings().collect()
    df_group_unmapped = group_selection.unmatched_ids().collect()
    df_group_edges = group_selection.edges().collect()

    assert df_group_mapping.columns == [
        "group_id",
        "input_id",
        "string_protein_id",
        "alias",
        "source",
    ]
    assert df_group_unmapped.columns == ["group_id", "input_id"]
    assert df_group_edges.columns == [
        "group_id",
        "string_id_a",
        "string_id_b",
        "combined_score",
    ]

    assert df_group_mapping.to_dicts() == [
        {
            "group_id": "G1",
            "input_id": "EGFR",
            "string_protein_id": "9606.ENSP0002",
            "alias": "EGFR",
            "source": "UniProt_GN_Name",
        },
        {
            "group_id": "G1",
            "input_id": "P04637",
            "string_protein_id": "9606.ENSP0001",
            "alias": "P04637",
            "source": "UniProt_AC",
        },
        {
            "group_id": "G1",
            "input_id": "P04637",
            "string_protein_id": "9606.ENSP9999",
            "alias": "P04637",
            "source": "UniProt_GN_Synonyms",
        },
        {
            "group_id": "G2",
            "input_id": "CDK2",
            "string_protein_id": "9606.ENSP0003",
            "alias": "CDK2",
            "source": "UniProt_GN_Name",
        },
        {
            "group_id": "G2",
            "input_id": "P04637",
            "string_protein_id": "9606.ENSP0001",
            "alias": "P04637",
            "source": "UniProt_AC",
        },
        {
            "group_id": "G2",
            "input_id": "P04637",
            "string_protein_id": "9606.ENSP9999",
            "alias": "P04637",
            "source": "UniProt_GN_Synonyms",
        },
    ]
    assert df_group_unmapped.to_dicts() == [
        {"group_id": "G1", "input_id": "MISSING"},
    ]
    assert df_group_edges.to_dicts() == [
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0002",
            "string_id_b": "9606.ENSP9999",
            "combined_score": 800,
        },
        {
            "group_id": "G2",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "group_id": "G2",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0003",
            "combined_score": 450,
        },
    ]


def test_group_selection_matches_equivalent_single_query_results(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)
    db = STRINGDatabase.from_files(aliases=file_aliases, links=file_links)

    group_selection = db.select_groups(
        {
            "G1": ["P04637", "EGFR", "MISSING"],
            "G2": ["P04637", "CDK2"],
        }
    ).with_min_combined_score(300)

    df_group_mapping = group_selection.mappings().collect()
    df_group_unmapped = group_selection.unmatched_ids().collect()
    df_group_edges = group_selection.edges().collect()
    selection_g1 = db.select_ids(["P04637", "EGFR", "MISSING"]).with_min_combined_score(
        300
    )
    selection_g2 = db.select_ids(["P04637", "CDK2"]).with_min_combined_score(300)

    df_single_mapping = pl.concat(
        [
            selection_g1.mappings()
            .collect()
            .with_columns(pl.lit("G1").alias("group_id"))
            .select(["group_id", "input_id", "string_protein_id", "alias", "source"]),
            selection_g2.mappings()
            .collect()
            .with_columns(pl.lit("G2").alias("group_id"))
            .select(["group_id", "input_id", "string_protein_id", "alias", "source"]),
        ]
    ).sort(["group_id", "input_id", "string_protein_id", "alias", "source"])
    df_single_unmapped = pl.concat(
        [
            selection_g1.unmatched_ids()
            .collect()
            .with_columns(pl.lit("G1").alias("group_id"))
            .select(["group_id", "input_id"]),
            selection_g2.unmatched_ids()
            .collect()
            .with_columns(pl.lit("G2").alias("group_id"))
            .select(["group_id", "input_id"]),
        ],
        how="vertical_relaxed",
    ).sort(["group_id", "input_id"])
    df_single_edges = pl.concat(
        [
            selection_g1.edges()
            .collect()
            .with_columns(pl.lit("G1").alias("group_id"))
            .select(["group_id", "string_id_a", "string_id_b", "combined_score"]),
            selection_g2.edges()
            .collect()
            .with_columns(pl.lit("G2").alias("group_id"))
            .select(["group_id", "string_id_a", "string_id_b", "combined_score"]),
        ]
    ).sort(["group_id", "string_id_a", "string_id_b"])

    assert df_group_mapping.filter(pl.col("group_id").is_in(["G1", "G2"])).equals(
        df_single_mapping
    )
    assert df_group_unmapped.equals(df_single_unmapped)
    assert df_group_edges.equals(df_single_edges)


def test_group_selection_reuses_cached_frames_and_recomputes_edges(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection_low = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_groups({"G1": ["P04637", "EGFR"], "G2": ["P04637", "CDK2"]})
        .with_min_combined_score(300)
    )

    df_mapping_low = selection_low.mappings().collect()
    df_unmapped_low = selection_low.unmatched_ids().collect()
    df_edges_low = selection_low.edges().collect()

    selection_high = selection_low.with_min_combined_score(600)

    assert selection_high.mappings().collect().equals(df_mapping_low)
    assert selection_high.unmatched_ids().collect().equals(df_unmapped_low)
    assert not selection_high.edges().collect().equals(df_edges_low)
    assert selection_high.edges().collect().to_dicts() == [
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 700,
        },
        {
            "group_id": "G1",
            "string_id_a": "9606.ENSP0002",
            "string_id_b": "9606.ENSP9999",
            "combined_score": 800,
        },
        {
            "group_id": "G2",
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0001",
            "combined_score": 999,
        },
    ]


def test_group_selection_resolves_normalized_ids_once_then_expands_membership(
    tmp_path: Path,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)

    selection = STRINGDatabase.from_files(
        aliases=file_aliases,
        links=file_links,
    ).select_groups(
        {
            "G1": ["sp|P04637|P53_HUMAN", "EGFR", "MISSING"],
            "G2": ["P04637", "CDK2", "MISSING"],
        }
    )

    assert selection._df_input_ids.to_dicts() == [  # pyright: ignore[reportPrivateUsage]
        {"input_id": "CDK2"},
        {"input_id": "EGFR"},
        {"input_id": "MISSING"},
        {"input_id": "P04637"},
    ]
    group_membership = selection._df_group_membership  # pyright: ignore[reportPrivateUsage]
    assert group_membership is not None
    assert group_membership.filter(
        pl.col("input_id").is_in(["MISSING", "P04637"])
    ).to_dicts() == [
        {"group_id": "G1", "input_id": "MISSING"},
        {"group_id": "G1", "input_id": "P04637"},
        {"group_id": "G2", "input_id": "MISSING"},
        {"group_id": "G2", "input_id": "P04637"},
    ]

    shared_mapping = (
        selection.mappings().collect().filter(pl.col("input_id") == "P04637")
    )
    assert shared_mapping.group_by("group_id").len().sort("group_id").to_dicts() == [
        {"group_id": "G1", "len": 2},
        {"group_id": "G2", "len": 2},
    ]
    assert selection.unmatched_ids().collect().to_dicts() == [
        {"group_id": "G1", "input_id": "MISSING"},
        {"group_id": "G2", "input_id": "MISSING"},
    ]
    assert selection.edges().collect().select("group_id", "string_id_b").to_dicts() == [
        {"group_id": "G1", "string_id_b": "9606.ENSP0001"},
        {"group_id": "G1", "string_id_b": "9606.ENSP0002"},
        {"group_id": "G1", "string_id_b": "9606.ENSP9999"},
        {"group_id": "G2", "string_id_b": "9606.ENSP0001"},
        {"group_id": "G2", "string_id_b": "9606.ENSP0003"},
    ]


def test_group_selection_builds_one_data_scan_per_cached_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(aliases=file_aliases, links=file_links)
    database = STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
    _ = database._alias_schema  # pyright: ignore[reportPrivateUsage]

    original_scan_aliases = stringdb_module.scan_aliases
    original_scan_links = stringdb_module.scan_links
    scan_counts = {"aliases": 0, "links": 0}

    def counted_scan_aliases(file_aliases: Path, version: str) -> pl.LazyFrame:
        scan_counts["aliases"] += 1
        return original_scan_aliases(file_aliases, version)

    def counted_scan_links(file_links: Path, version: str) -> pl.LazyFrame:
        scan_counts["links"] += 1
        return original_scan_links(file_links, version)

    monkeypatch.setattr(stringdb_module, "scan_aliases", counted_scan_aliases)
    monkeypatch.setattr(stringdb_module, "scan_links", counted_scan_links)

    selection = database.select_groups(
        {
            "G1": ["sp|P04637|P53_HUMAN", "EGFR"],
            "G2": ["P04637", "CDK2"],
        }
    )
    assert selection.mappings().collect().equals(selection.mappings().collect())
    assert (
        selection.unmatched_ids().collect().equals(selection.unmatched_ids().collect())
    )
    assert selection.edges().collect().equals(selection.edges().collect())
    # Taxon compatibility validation adds one bounded source scan, but the
    # count remains independent of the number of groups.
    assert all(value >= 1 for value in scan_counts.values())


def test_direct_string_namespace_uses_links_without_aliases(tmp_path: Path) -> None:
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n9606.ENSP0001 9606.ENSP0002 500\n",
        should_gzip=False,
    )

    selection = STRINGDatabase.from_files(links=file_links).select_ids(
        ["9606.ENSP0001", "9606.ENSP0002"], namespace="string"
    )

    assert selection.edges().collect().to_dicts() == [
        {
            "string_id_a": "9606.ENSP0001",
            "string_id_b": "9606.ENSP0002",
            "combined_score": 500,
        }
    ]
    assert selection.unmatched_ids().collect().to_dicts() == []
    with pytest.raises(
        (CapabilityError, pl.exceptions.ComputeError), match="namespace='alias'"
    ):
        selection.mappings().collect()


def test_canonical_edge_score_conflict_fails_before_threshold(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_text_or_gzip(
        file_aliases,
        "string_protein_id\talias\tsource\n"
        "9606.ENSP0001\tTP53\tUniProt_GN_Name\n"
        "9606.ENSP0002\tEGFR\tUniProt_GN_Name\n",
        should_gzip=False,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n"
        "9606.ENSP0001 9606.ENSP0002 400\n"
        "9606.ENSP0002 9606.ENSP0001 700\n",
        should_gzip=False,
    )

    selection = (
        STRINGDatabase.from_files(aliases=file_aliases, links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_min_combined_score(700)
    )
    with pytest.raises(
        (IntegrityError, pl.exceptions.ComputeError),
        match="conflicting combined_score",
    ):
        selection.edges().collect()


def test_taxon_mismatch_is_rejected_and_gzip_inputs_warn(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    file_aliases = tmp_path / "aliases.txt.gz"
    file_links = tmp_path / "links.txt.gz"
    _write_text_or_gzip(
        file_aliases,
        "string_protein_id\talias\tsource\n9606.ENSP0001\tTP53\tUniProt_GN_Name\n",
        should_gzip=True,
    )
    _write_text_or_gzip(
        file_links,
        "protein1 protein2 combined_score\n10090.ENSP0001 10090.ENSP0002 500\n",
        should_gzip=True,
    )

    caplog.set_level("WARNING")
    selection = STRINGDatabase.from_files(
        aliases=file_aliases, links=file_links
    ).select_ids(["TP53"])
    assert "gzip-compressed" in caplog.text
    with pytest.raises(
        (IntegrityError, pl.exceptions.ComputeError), match="incompatible taxon"
    ):
        selection.edges().collect()
