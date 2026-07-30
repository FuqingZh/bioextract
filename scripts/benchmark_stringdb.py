from __future__ import annotations

import argparse
import json
import platform
import random
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

import polars as pl

from bioextract.stringdb import STRINGDatabase

try:
    import resource
except ImportError:  # pragma: no cover
    resource = None


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="benchmark_stringdb",
        description=(
            "Benchmark repeated single-query scans versus grouped STRING edge "
            "extraction on synthetic data."
        ),
    )
    parser.add_argument("--num_proteins", type=int, default=20_000)
    parser.add_argument("--num_links", type=int, default=400_000)
    parser.add_argument("--query_size", type=int, default=500)
    parser.add_argument("--num_groups", type=int, default=8)
    parser.add_argument("--thr_score_min", type=int, default=400)
    parser.add_argument("--seed", type=int, default=17)
    return parser


def _write_aliases(file_aliases: Path, *, num_proteins: int) -> list[str]:
    protein_ids: list[str] = []
    with file_aliases.open("w", encoding="utf-8") as handle:
        handle.write("#string_protein_id\talias\tsource\n")
        for idx in range(num_proteins):
            string_id = f"9606.ENSP{idx:011d}"
            protein_id = f"P{idx:05d}"
            protein_ids.append(protein_id)
            handle.write(f"{string_id}\t{protein_id}\tUniProt_AC\n")
            handle.write(f"{string_id}\tGENE{idx:05d}\tUniProt_GN_Name\n")
    return protein_ids


def _write_links(
    file_links: Path,
    *,
    num_proteins: int,
    num_links: int,
    seed: int,
) -> None:
    rng = random.Random(seed)
    with file_links.open("w", encoding="utf-8") as handle:
        handle.write("protein1 protein2 combined_score\n")
        for _ in range(num_links):
            idx_a = rng.randrange(num_proteins)
            idx_b = rng.randrange(num_proteins)
            score = rng.randint(150, 999)
            handle.write(f"9606.ENSP{idx_a:011d} 9606.ENSP{idx_b:011d} {score}\n")


def _create_group_queries(
    *,
    protein_ids: list[str],
    num_groups: int,
    query_size: int,
) -> dict[str, list[str]]:
    if num_groups <= 0:
        return {}

    num_proteins = len(protein_ids)
    ids_by_group: dict[str, list[str]] = {}
    for idx_group in range(num_groups):
        ids_group = [
            protein_ids[(idx_group * query_size + idx_query) % num_proteins]
            for idx_query in range(query_size)
        ]
        ids_by_group[f"G{idx_group + 1:03d}"] = ids_group
    return ids_by_group


def _extract_edges_repeated_single_queries(
    *,
    db: STRINGDatabase,
    ids_by_group: dict[str, list[str]],
    thr_score_min: int,
) -> pl.DataFrame:
    df_edges_grouped: list[pl.DataFrame] = []
    for group_id, input_ids in ids_by_group.items():
        df_edges_group = (
            db.select_ids(input_ids)
            .with_min_combined_score(thr_score_min)
            .extract_edges()
            .with_columns(pl.lit(group_id).alias("GroupId"))
            .select(["GroupId", "StringIdA", "StringIdB", "Score"])
        )
        df_edges_grouped.append(df_edges_group)

    if not df_edges_grouped:
        return pl.DataFrame(
            schema={
                "GroupId": pl.String,
                "StringIdA": pl.String,
                "StringIdB": pl.String,
                "Score": pl.Int64,
            }
        )

    return pl.concat(df_edges_grouped).sort(["GroupId", "StringIdA", "StringIdB"])


def _run_case(label: str, fn_case: Callable[[], pl.DataFrame]) -> dict[str, object]:
    time_started = time.perf_counter()
    df_edges = fn_case()
    elapsed_seconds = time.perf_counter() - time_started
    peak_rss_mb = None
    if resource is not None:
        peak_rss_raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_rss_mb = peak_rss_raw / 1024.0

    return {
        "label": label,
        "wall_time_s": round(elapsed_seconds, 4),
        "num_edges": int(df_edges.height),
        "peak_rss_mb_process_global": (
            None if peak_rss_mb is None else round(peak_rss_mb, 2)
        ),
    }


def main() -> None:
    args = _build_parser().parse_args()

    with tempfile.TemporaryDirectory(prefix="bioextract-bench-") as dir_tmp:
        dir_tmp_path = Path(dir_tmp)
        file_aliases = dir_tmp_path / "aliases.txt"
        file_links = dir_tmp_path / "links.txt"

        protein_ids = _write_aliases(file_aliases, num_proteins=args.num_proteins)
        _write_links(
            file_links,
            num_proteins=args.num_proteins,
            num_links=args.num_links,
            seed=args.seed,
        )
        ids_by_group = _create_group_queries(
            protein_ids=protein_ids,
            num_groups=args.num_groups,
            query_size=args.query_size,
        )
        db = STRINGDatabase.from_files(
            aliases=file_aliases,
            links=file_links,
        )

        case_repeated = _run_case(
            "repeated_single_queries",
            lambda: _extract_edges_repeated_single_queries(
                db=db,
                ids_by_group=ids_by_group,
                thr_score_min=args.thr_score_min,
            ),
        )
        case_grouped = _run_case(
            "grouped_query",
            lambda: (
                db.select_groups(ids_by_group)
                .with_min_combined_score(args.thr_score_min)
                .extract_edges()
            ),
        )

        report = {
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "polars": pl.__version__,
            },
            "inputs": {
                "num_proteins": args.num_proteins,
                "num_links": args.num_links,
                "num_groups": args.num_groups,
                "query_size": args.query_size,
                "total_input_ids": args.num_groups * args.query_size,
                "thr_score_min": args.thr_score_min,
                "seed": args.seed,
                "file_aliases_bytes": file_aliases.stat().st_size,
                "file_links_bytes": file_links.stat().st_size,
                "rows_scanned_aliases": args.num_proteins * 2,
                "rows_scanned_links": args.num_links,
            },
            "results": [case_repeated, case_grouped],
            "notes": [
                "peak_rss_mb_process_global is process-global and not isolated per case",
                "grouped_query scans aliases and links once for all groups",
                "repeated_single_queries reruns single-query extraction per group",
            ],
        }
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
