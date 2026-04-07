from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from bioextract.stringdb import StringDb, StringResourceLimits


def _write_text_or_gzip(file_out: Path, content: str, *, should_gzip: bool) -> None:
    if should_gzip:
        with gzip.open(file_out, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return

    file_out.write_text(content, encoding="utf-8")


def _write_demo_string_files(
    *,
    file_aliases: Path,
    file_links: Path,
    should_gzip: bool = False,
) -> None:
    _write_text_or_gzip(
        file_aliases,
        "\n".join(
            [
                "#string_protein_id\talias\tsource",
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
        file_links,
        "\n".join(
            [
                "protein1 protein2 combined_score",
                "9606.ENSP0001 9606.ENSP0002 400",
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

    selection = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
    ).select_ids([" sp|P04637|P53_HUMAN ", " EGFR ", "  "])

    df_result = selection.extract_string_mapping()

    assert df_result.columns == ["InputId", "StringId", "MapSource"]
    assert df_result.to_dicts() == [
        {
            "InputId": "EGFR",
            "StringId": "9606.ENSP0002",
            "MapSource": "UniProt_GN_Name",
        },
        {
            "InputId": "P04637",
            "StringId": "9606.ENSP0001",
            "MapSource": "UniProt_AC",
        },
    ]


def test_stringdb_single_query_smoke(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    result = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_score_min(400)
        .extract_result()
    )

    assert result.df_edges.columns == ["StringIdA", "StringIdB", "Score"]
    assert result.df_protein_map.columns == ["InputId", "StringId", "MapSource"]
    assert result.df_unmapped.columns == ["InputId"]
    assert result.metrics.num_input_ids == 2
    assert result.metrics.num_edges == 1


def test_stringdb_group_query_smoke(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    selection = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_groups({"G1": ["TP53", "EGFR"], "G2": ["P04637", "CDK2"]})
        .with_score_min(400)
    )

    df_edges = selection.extract_edges()
    df_metrics = selection.extract_metrics()

    assert df_edges.columns == ["GroupId", "StringIdA", "StringIdB", "Score"]
    assert df_metrics.columns == [
        "GroupId",
        "NumInputIds",
        "NumMappedIds",
        "NumUnmappedIds",
        "NumStringIds",
        "NumEdges",
    ]
    assert df_metrics.to_dicts() == [
        {
            "GroupId": "G1",
            "NumInputIds": 2,
            "NumMappedIds": 2,
            "NumUnmappedIds": 0,
            "NumStringIds": 2,
            "NumEdges": 1,
        },
        {
            "GroupId": "G2",
            "NumInputIds": 2,
            "NumMappedIds": 2,
            "NumUnmappedIds": 0,
            "NumStringIds": 3,
            "NumEdges": 1,
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

        selection = StringDb.from_files(
            file_aliases=file_aliases,
            file_links=file_links,
        ).select_ids(["TP53", "EGFR"])

        df_result = selection.extract_string_mapping()
        assert df_result.to_dicts() == [
            {
                "InputId": "EGFR",
                "StringId": "9606.ENSP0002",
                "MapSource": "UniProt_GN_Name",
            },
            {
                "InputId": "TP53",
                "StringId": "9606.ENSP0001",
                "MapSource": "UniProt_GN_Name",
            },
        ]


def test_extract_result_returns_bundle_and_metrics(tmp_path: Path) -> None:
    for should_gzip in (False, True):
        suffix = ".txt.gz" if should_gzip else ".txt"
        file_aliases = tmp_path / f"aliases{suffix}"
        file_links = tmp_path / f"links{suffix}"
        _write_demo_string_files(
            file_aliases=file_aliases,
            file_links=file_links,
            should_gzip=should_gzip,
        )

        result = (
            StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
            .select_ids(["sp|P04637|P53_HUMAN", "EGFR", "CDK2", "MISSING"])
            .with_score_min(300)
            .extract_result()
        )

        assert result.df_edges.columns == ["StringIdA", "StringIdB", "Score"]
        assert result.df_protein_map.columns == ["InputId", "StringId", "MapSource"]
        assert result.df_unmapped.columns == ["InputId"]

        assert result.df_edges.to_dicts() == [
            {
                "StringIdA": "9606.ENSP0001",
                "StringIdB": "9606.ENSP0002",
                "Score": 700,
            },
            {
                "StringIdA": "9606.ENSP0001",
                "StringIdB": "9606.ENSP0003",
                "Score": 450,
            },
            {
                "StringIdA": "9606.ENSP0002",
                "StringIdB": "9606.ENSP9999",
                "Score": 800,
            },
        ]
        assert result.df_protein_map.to_dicts() == [
            {
                "InputId": "CDK2",
                "StringId": "9606.ENSP0003",
                "MapSource": "UniProt_GN_Name",
            },
            {
                "InputId": "EGFR",
                "StringId": "9606.ENSP0002",
                "MapSource": "UniProt_GN_Name",
            },
            {
                "InputId": "P04637",
                "StringId": "9606.ENSP0001",
                "MapSource": "UniProt_AC",
            },
            {
                "InputId": "P04637",
                "StringId": "9606.ENSP9999",
                "MapSource": "UniProt_GN_Synonyms",
            },
        ]
        assert result.df_unmapped.to_dicts() == [{"InputId": "MISSING"}]

        assert result.metrics.num_input_ids == 4
        assert result.metrics.num_mapped_ids == 3
        assert result.metrics.num_unmapped_ids == 1
        assert result.metrics.num_string_ids == 4
        assert result.metrics.num_edges == 3


def test_extract_result_handles_empty_input_set(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    result = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_ids([])
        .extract_result()
    )

    assert result.df_edges.schema == {
        "StringIdA": pl.String,
        "StringIdB": pl.String,
        "Score": pl.Int64,
    }
    assert result.df_protein_map.schema == {
        "InputId": pl.String,
        "StringId": pl.String,
        "MapSource": pl.String,
    }
    assert result.df_unmapped.schema == {"InputId": pl.String}
    assert result.metrics.num_input_ids == 0
    assert result.metrics.num_mapped_ids == 0
    assert result.metrics.num_unmapped_ids == 0
    assert result.metrics.num_string_ids == 0
    assert result.metrics.num_edges == 0


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
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_ids(["TP53"])
        .extract_string_mapping()
    )
    df_custom = (
        StringDb.from_files(
            file_aliases=file_aliases,
            file_links=file_links,
            source_rank_map={"CustomSource": 1, **StringDb.DEFAULT_SOURCE_RANK_MAP},
        )
        .select_ids(["TP53"])
        .extract_string_mapping()
    )

    assert df_default.to_dicts() == [
        {
            "InputId": "TP53",
            "StringId": "9606.ENSP0001",
            "MapSource": "UniProt_GN_Name",
        }
    ]
    assert df_custom.to_dicts() == [
        {
            "InputId": "TP53",
            "StringId": "9606.ENSP0001",
            "MapSource": "CustomSource",
        }
    ]


def test_selection_reuses_cached_frames_and_result(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    selection = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_score_min(400)
    )

    df_protein_map_first = selection.extract_string_mapping()
    df_protein_map_second = selection.extract_string_mapping()
    assert df_protein_map_first is df_protein_map_second

    df_edges_first = selection.extract_edges()
    df_edges_second = selection.extract_edges()
    assert df_edges_first is df_edges_second

    result_first = selection.extract_result()
    result_second = selection.extract_result()
    assert result_first is result_second
    assert result_first.df_protein_map is df_protein_map_first
    assert result_first.df_edges is df_edges_first


def test_with_score_min_reuses_map_cache_but_recomputes_edges(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    selection_low = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_ids(["TP53", "EGFR"])
        .with_score_min(300)
    )
    selection_low_map = selection_low.extract_string_mapping()
    selection_low_unmapped = selection_low.extract_unmapped_input_ids()
    selection_low_edges = selection_low.extract_edges()

    selection_high = selection_low.with_score_min(600)

    assert selection_high.extract_string_mapping() is selection_low_map
    assert selection_high.extract_unmapped_input_ids() is selection_low_unmapped
    assert selection_high.extract_edges() is not selection_low_edges
    assert selection_low.extract_edges().to_dicts() == [
        {
            "StringIdA": "9606.ENSP0001",
            "StringIdB": "9606.ENSP0002",
            "Score": 700,
        }
    ]
    assert selection_high.extract_edges().to_dicts() == [
        {
            "StringIdA": "9606.ENSP0001",
            "StringIdB": "9606.ENSP0002",
            "Score": 700,
        }
    ]


def test_from_files_rejects_files_over_limit(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    with pytest.raises(ValueError, match="aliases file"):
        StringDb.from_files(
            file_aliases=file_aliases,
            file_links=file_links,
            limits=StringResourceLimits(file_aliases_bytes_max=1),
        )

    with pytest.raises(ValueError, match="links file"):
        StringDb.from_files(
            file_aliases=file_aliases,
            file_links=file_links,
            limits=StringResourceLimits(file_links_bytes_max=1),
        )


def test_from_files_rejects_missing_files(tmp_path: Path) -> None:
    file_aliases = tmp_path / "missing_aliases.txt"
    file_links = tmp_path / "missing_links.txt"

    with pytest.raises(FileNotFoundError, match="aliases file not found"):
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)

    file_aliases.write_text("", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="links file not found"):
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)


def test_from_files_rejects_unsupported_version(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    with pytest.raises(ValueError, match="Unsupported STRING version"):
        StringDb.from_files(
            file_aliases=file_aliases,
            file_links=file_links,
            version="v11.5",  # type: ignore[arg-type]
        )


def test_select_ids_rejects_input_count_over_limit(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    db = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
        limits=StringResourceLimits(num_input_ids_max=1),
    )

    with pytest.raises(ValueError, match="Normalized input ID count exceeds"):
        db.select_ids(["TP53", "EGFR"])


def test_select_groups_rejects_invalid_group_shape(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    db = StringDb.from_files(file_aliases=file_aliases, file_links=file_links)

    with pytest.raises(ValueError, match="non-empty"):
        db.select_groups({"  ": ["TP53"]})

    with pytest.raises(ValueError, match="unique after normalization"):
        db.select_groups({"A": ["TP53"], " A ": ["EGFR"]})


def test_select_groups_rejects_group_and_input_limits(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    db_group_limit = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
        limits=StringResourceLimits(num_groups_max=1),
    )
    with pytest.raises(ValueError, match="Group count exceeds"):
        db_group_limit.select_groups({"G1": ["TP53"], "G2": ["EGFR"]})

    db_input_limit = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
        limits=StringResourceLimits(num_input_ids_max=2),
    )
    with pytest.raises(ValueError, match="Normalized input ID count exceeds"):
        db_input_limit.select_groups({"G1": ["TP53", "EGFR"], "G2": ["CDK2"]})


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

    selection = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
    ).select_ids(["TP53"])

    with pytest.raises(ValueError, match="aliases file"):
        selection.extract_string_mapping()


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

    selection = StringDb.from_files(
        file_aliases=file_aliases,
        file_links=file_links,
    ).select_ids(["TP53"])

    with pytest.raises(ValueError, match="links file"):
        selection.extract_edges()


def test_group_selection_extracts_flat_outputs_and_metrics(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    group_selection = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_groups(
            {
                "G1": ["sp|P04637|P53_HUMAN", "EGFR", "MISSING"],
                "G2": ["P04637", "CDK2"],
                "G3": [],
            }
        )
        .with_score_min(300)
    )

    df_group_mapping = group_selection.extract_string_mapping()
    df_group_unmapped = group_selection.extract_unmapped_input_ids()
    df_group_edges = group_selection.extract_edges()
    df_group_metrics = group_selection.extract_metrics()

    assert df_group_mapping.columns == ["GroupId", "InputId", "StringId", "MapSource"]
    assert df_group_unmapped.columns == ["GroupId", "InputId"]
    assert df_group_edges.columns == ["GroupId", "StringIdA", "StringIdB", "Score"]
    assert df_group_metrics.columns == [
        "GroupId",
        "NumInputIds",
        "NumMappedIds",
        "NumUnmappedIds",
        "NumStringIds",
        "NumEdges",
    ]

    assert df_group_mapping.to_dicts() == [
        {
            "GroupId": "G1",
            "InputId": "EGFR",
            "StringId": "9606.ENSP0002",
            "MapSource": "UniProt_GN_Name",
        },
        {
            "GroupId": "G1",
            "InputId": "P04637",
            "StringId": "9606.ENSP0001",
            "MapSource": "UniProt_AC",
        },
        {
            "GroupId": "G1",
            "InputId": "P04637",
            "StringId": "9606.ENSP9999",
            "MapSource": "UniProt_GN_Synonyms",
        },
        {
            "GroupId": "G2",
            "InputId": "CDK2",
            "StringId": "9606.ENSP0003",
            "MapSource": "UniProt_GN_Name",
        },
        {
            "GroupId": "G2",
            "InputId": "P04637",
            "StringId": "9606.ENSP0001",
            "MapSource": "UniProt_AC",
        },
        {
            "GroupId": "G2",
            "InputId": "P04637",
            "StringId": "9606.ENSP9999",
            "MapSource": "UniProt_GN_Synonyms",
        },
    ]
    assert df_group_unmapped.to_dicts() == [
        {"GroupId": "G1", "InputId": "MISSING"},
    ]
    assert df_group_edges.to_dicts() == [
        {
            "GroupId": "G1",
            "StringIdA": "9606.ENSP0001",
            "StringIdB": "9606.ENSP0002",
            "Score": 700,
        },
        {
            "GroupId": "G1",
            "StringIdA": "9606.ENSP0002",
            "StringIdB": "9606.ENSP9999",
            "Score": 800,
        },
        {
            "GroupId": "G2",
            "StringIdA": "9606.ENSP0001",
            "StringIdB": "9606.ENSP0003",
            "Score": 450,
        },
    ]
    assert df_group_metrics.to_dicts() == [
        {
            "GroupId": "G1",
            "NumInputIds": 3,
            "NumMappedIds": 2,
            "NumUnmappedIds": 1,
            "NumStringIds": 3,
            "NumEdges": 2,
        },
        {
            "GroupId": "G2",
            "NumInputIds": 2,
            "NumMappedIds": 2,
            "NumUnmappedIds": 0,
            "NumStringIds": 3,
            "NumEdges": 1,
        },
        {
            "GroupId": "G3",
            "NumInputIds": 0,
            "NumMappedIds": 0,
            "NumUnmappedIds": 0,
            "NumStringIds": 0,
            "NumEdges": 0,
        },
    ]


def test_group_selection_matches_equivalent_single_query_results(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)
    db = StringDb.from_files(file_aliases=file_aliases, file_links=file_links)

    group_selection = db.select_groups(
        {
            "G1": ["P04637", "EGFR", "MISSING"],
            "G2": ["P04637", "CDK2"],
        }
    ).with_score_min(300)

    df_group_mapping = group_selection.extract_string_mapping()
    df_group_unmapped = group_selection.extract_unmapped_input_ids()
    df_group_edges = group_selection.extract_edges()
    df_group_metrics = group_selection.extract_metrics()

    result_g1 = db.select_ids(["P04637", "EGFR", "MISSING"]).with_score_min(
        300
    ).extract_result()
    result_g2 = db.select_ids(["P04637", "CDK2"]).with_score_min(300).extract_result()

    df_single_mapping = pl.concat(
        [
            result_g1.df_protein_map.with_columns(pl.lit("G1").alias("GroupId")).select(
                ["GroupId", "InputId", "StringId", "MapSource"]
            ),
            result_g2.df_protein_map.with_columns(pl.lit("G2").alias("GroupId")).select(
                ["GroupId", "InputId", "StringId", "MapSource"]
            ),
        ]
    ).sort(["GroupId", "InputId", "StringId", "MapSource"])
    df_single_unmapped = pl.concat(
        [
            result_g1.df_unmapped.with_columns(pl.lit("G1").alias("GroupId")).select(
                ["GroupId", "InputId"]
            ),
            result_g2.df_unmapped.with_columns(pl.lit("G2").alias("GroupId")).select(
                ["GroupId", "InputId"]
            ),
        ],
        how="vertical_relaxed",
    ).sort(["GroupId", "InputId"])
    df_single_edges = pl.concat(
        [
            result_g1.df_edges.with_columns(pl.lit("G1").alias("GroupId")).select(
                ["GroupId", "StringIdA", "StringIdB", "Score"]
            ),
            result_g2.df_edges.with_columns(pl.lit("G2").alias("GroupId")).select(
                ["GroupId", "StringIdA", "StringIdB", "Score"]
            ),
        ]
    ).sort(["GroupId", "StringIdA", "StringIdB"])
    df_single_metrics = pl.DataFrame(
        [
            {
                "GroupId": "G1",
                "NumInputIds": result_g1.metrics.num_input_ids,
                "NumMappedIds": result_g1.metrics.num_mapped_ids,
                "NumUnmappedIds": result_g1.metrics.num_unmapped_ids,
                "NumStringIds": result_g1.metrics.num_string_ids,
                "NumEdges": result_g1.metrics.num_edges,
            },
            {
                "GroupId": "G2",
                "NumInputIds": result_g2.metrics.num_input_ids,
                "NumMappedIds": result_g2.metrics.num_mapped_ids,
                "NumUnmappedIds": result_g2.metrics.num_unmapped_ids,
                "NumStringIds": result_g2.metrics.num_string_ids,
                "NumEdges": result_g2.metrics.num_edges,
            },
        ],
        schema={
            "GroupId": pl.String,
            "NumInputIds": pl.Int64,
            "NumMappedIds": pl.Int64,
            "NumUnmappedIds": pl.Int64,
            "NumStringIds": pl.Int64,
            "NumEdges": pl.Int64,
        },
    ).sort("GroupId")

    assert df_group_mapping.filter(pl.col("GroupId").is_in(["G1", "G2"])).equals(
        df_single_mapping
    )
    assert df_group_unmapped.equals(df_single_unmapped)
    assert df_group_edges.equals(df_single_edges)
    assert df_group_metrics.filter(pl.col("GroupId").is_in(["G1", "G2"])).equals(
        df_single_metrics
    )


def test_group_selection_reuses_cached_frames_and_recomputes_edges(tmp_path: Path) -> None:
    file_aliases = tmp_path / "aliases.txt"
    file_links = tmp_path / "links.txt"
    _write_demo_string_files(file_aliases=file_aliases, file_links=file_links)

    selection_low = (
        StringDb.from_files(file_aliases=file_aliases, file_links=file_links)
        .select_groups({"G1": ["P04637", "EGFR"], "G2": ["P04637", "CDK2"]})
        .with_score_min(300)
    )

    df_mapping_low = selection_low.extract_string_mapping()
    df_unmapped_low = selection_low.extract_unmapped_input_ids()
    df_edges_low = selection_low.extract_edges()

    selection_high = selection_low.with_score_min(600)

    assert selection_high.extract_string_mapping() is df_mapping_low
    assert selection_high.extract_unmapped_input_ids() is df_unmapped_low
    assert selection_high.extract_edges() is not df_edges_low
    assert selection_high.extract_metrics().to_dicts() == [
        {
            "GroupId": "G1",
            "NumInputIds": 2,
            "NumMappedIds": 2,
            "NumUnmappedIds": 0,
            "NumStringIds": 3,
            "NumEdges": 2,
        },
        {
            "GroupId": "G2",
            "NumInputIds": 2,
            "NumMappedIds": 2,
            "NumUnmappedIds": 0,
            "NumStringIds": 3,
            "NumEdges": 0,
        },
    ]
