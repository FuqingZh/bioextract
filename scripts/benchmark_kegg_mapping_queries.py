"""Benchmark public nested KEGG selections against a temporary flat lookup."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from bioextract import KEGGDatabase
from bioextract.kegg.mapping.constant import KEGGNamespace
from bioextract.kegg.mapping.query import KeggSelection


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("publication", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--organism", action="append", default=[])
    args = parser.parse_args()

    database = KEGGDatabase.from_duckdb(args.publication)
    if args.organism:
        database = database.with_organisms(args.organism)
    records: list[dict[str, Any]] = []
    with duckdb.connect() as connection:
        publication_literal = "'" + str(args.publication).replace("'", "''") + "'"
        connection.execute(f"ATTACH {publication_literal} AS publication (READ_ONLY)")
        namespaces: dict[KEGGNamespace, tuple[str, str]] = {
            "uniprot": (
                "uniprot_mappings",
                "uniprot_id",
            ),
            "ncbi_gene": (
                "ncbi_gene_mappings",
                "ncbi_gene_id",
            ),
        }
        for namespace, (list_column, value_column) in namespaces.items():
            where = ""
            parameters: list[str] = []
            if args.organism:
                where = (
                    "WHERE organism_code IN ("
                    + ", ".join("?" for _ in args.organism)
                    + ")"
                )
                parameters.extend(args.organism)
            lookup_table = f"benchmark_lookup_{namespace}"
            lookup_started = time.perf_counter()
            connection.execute(
                f"CREATE TEMP TABLE {lookup_table} AS "
                f"SELECT DISTINCT mapping.{value_column} AS input_id, "
                "organism_code, kegg_gene_id "
                "FROM publication.gene_annotation, "
                f"UNNEST({list_column}) AS item(mapping) {where}",
                parameters,
            )
            lookup_build_seconds = time.perf_counter() - lookup_started
            available = [
                str(row[0])
                for row in connection.execute(
                    f"SELECT DISTINCT mapping.{value_column} "
                    "FROM publication.gene_annotation, "
                    f"UNNEST({list_column}) AS item(mapping) "
                    f"{where} ORDER BY mapping.{value_column}",
                    parameters,
                ).fetchall()
            ]
            for requested in (1, 100, 10_000):
                ids = available[:requested]
                if not ids:
                    continue
                selection: KeggSelection = database.select_ids(ids, namespace=namespace)
                input_frame = pl.DataFrame({"input_id": ids})
                connection.register("benchmark_input", input_frame)
                for run in ("cold", "warm"):
                    started = time.perf_counter()
                    frame = selection.matches().collect()
                    elapsed = time.perf_counter() - started
                    checksum = hashlib.sha256(
                        "\n".join(
                            "\t".join(str(value) for value in row)
                            for row in frame.sort(frame.columns).iter_rows()
                        ).encode()
                    ).hexdigest()
                    records.append(
                        {
                            "namespace": namespace,
                            "requested_input_count": requested,
                            "actual_input_count": len(ids),
                            "run": run,
                            "elapsed_seconds": elapsed,
                            "match_count": frame.height,
                            "checksum": checksum,
                            "implementation": "aggregate_public_api",
                        }
                    )
                    started = time.perf_counter()
                    baseline_rows = connection.execute(
                        f"SELECT input.input_id, '{namespace}', gene.organism_code, "
                        "gene.kegg_gene_id FROM benchmark_input AS input "
                        f"JOIN {lookup_table} AS gene USING (input_id) "
                        "ORDER BY ALL"
                    ).fetchall()
                    baseline_elapsed = time.perf_counter() - started
                    baseline_checksum = hashlib.sha256(
                        "\n".join(
                            "\t".join(str(value) for value in row)
                            for row in baseline_rows
                        ).encode()
                    ).hexdigest()
                    records.append(
                        {
                            "namespace": namespace,
                            "requested_input_count": requested,
                            "actual_input_count": len(ids),
                            "run": run,
                            "elapsed_seconds": baseline_elapsed,
                            "match_count": len(baseline_rows),
                            "checksum": baseline_checksum,
                            "implementation": "temporary_flat_sql",
                            "lookup_build_seconds": lookup_build_seconds,
                        }
                    )
                    if baseline_checksum != checksum:
                        raise RuntimeError(
                            f"KEGG selection checksum mismatch: {namespace}, {requested}"
                        )
                connection.unregister("benchmark_input")

    report = {
        "publication": str(args.publication),
        "organism_scope": sorted(set(args.organism)) or None,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "records": records,
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
