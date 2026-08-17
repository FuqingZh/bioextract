"""Benchmark representative public lazy-relation execution scenarios."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import statistics
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

import duckdb
import polars as pl

from bioextract import (
    ChEBIDatabase,
    KEGGDatabase,
    RheaDatabase,
    UniProtDatabase,
    inspect_publication,
)
from bioextract.uniprot._query import UniProtSelection

_CASE_NAMES = (
    "idmapping_narrow",
    "idmapping_wide",
    "rhea_matches",
    "rhea_bundle",
    "uniprot_proteins",
    "uniprot_bundle",
    "chebi_matches",
    "chebi_bundle",
    "kegg_matches",
    "kegg_reactions_narrow",
    "kegg_reactions_wide",
    "kegg_bundle",
)
_THREAD_ENVIRONMENT = (
    "BIOEXTRACT_TEST_THREADS",
    "POLARS_MAX_THREADS",
    "RAYON_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
)


@dataclass(frozen=True, slots=True)
class BenchmarkInputs:
    """Caller-supplied publications and identifiers for one benchmark run."""

    idmapping: Path
    uniprot_kb: Path
    rhea: Path
    chebi: Path
    kegg: Path
    taxon_id: str
    protein_ids: tuple[str, ...]
    kegg_compound_id: str


type CaseOperation = Callable[[], Mapping[str, pl.DataFrame]]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--idmapping", type=Path, required=True)
    parser.add_argument("--uniprot-kb", type=Path, required=True)
    parser.add_argument("--rhea", type=Path, required=True)
    parser.add_argument("--chebi", type=Path, required=True)
    parser.add_argument("--kegg", type=Path, required=True)
    parser.add_argument("--taxon-id", required=True)
    parser.add_argument("--protein-id", action="append", required=True)
    parser.add_argument("--kegg-compound-id", required=True)
    parser.add_argument("--case", action="append", choices=_CASE_NAMES)
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--label", default="working-tree")
    parser.add_argument("--output", type=Path)
    return parser


def _collect(frame: pl.LazyFrame, *columns: str) -> pl.DataFrame:
    return frame.select(*columns).collect(engine="streaming")


def _idmapping_narrow(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    database = UniProtDatabase.from_duckdb(inputs.idmapping)
    frame = (
        database.scan_mapping(taxon_ids=[inputs.taxon_id])
        .select("uniprot_id")
        .filter(pl.col("uniprot_id").is_in(inputs.protein_ids))
        .collect(engine="streaming")
    )
    return {"idmapping": frame}


def _idmapping_wide(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    database = UniProtDatabase.from_duckdb(inputs.idmapping)
    frame = (
        database.scan_mapping(taxon_ids=[inputs.taxon_id])
        .filter(pl.col("uniprot_id").is_in(inputs.protein_ids))
        .collect(engine="streaming")
    )
    return {"idmapping": frame}


def _rhea_matches(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = RheaDatabase.from_duckdb(inputs.rhea).select_reactions(
        inputs.protein_ids,
        namespace="uniprot",
    )
    return {"matches": _collect(selection.matches(), "input_id")}


def _rhea_bundle(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = RheaDatabase.from_duckdb(inputs.rhea).select_reactions(
        inputs.protein_ids,
        namespace="uniprot",
    )
    return {
        "matches": _collect(selection.matches(), "input_id"),
        "reactions": _collect(
            selection.reactions(),
            "input_id",
            "rhea_id",
            "master_id",
            "direction",
        ),
        "participants": _collect(
            selection.participants(),
            "input_id",
            "rhea_id",
            "master_id",
            "direction",
            "side",
            "directional_role",
            "chebi_id",
            "coefficient_numeric",
        ),
        "unmatched": _collect(
            selection.unmatched_ids(),
            "input_id",
            "input_namespace",
        ),
    }


def _uniprot_selection(inputs: BenchmarkInputs):
    return cast(
        UniProtSelection,
        UniProtDatabase.from_duckdb(inputs.uniprot_kb).select_ids(
            inputs.protein_ids,
            namespace="uniprot",
            taxon_ids=[inputs.taxon_id],
        ),
    )


def _uniprot_proteins(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _uniprot_selection(inputs)
    return {"proteins": _collect(selection.proteins(), "input_id")}


def _uniprot_bundle(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _uniprot_selection(inputs)
    return {
        "proteins": _collect(selection.proteins(), "input_id"),
        "ec_numbers": _collect(
            selection.ec_numbers(),
            "input_id",
            "ec_number",
        ),
        "unmatched": _collect(
            selection.unmatched_ids(),
            "input_id",
            "input_namespace",
            "reason",
        ),
    }


def _chebi_selection(inputs: BenchmarkInputs):
    return ChEBIDatabase.from_duckdb(inputs.chebi).select_compounds(
        [inputs.kegg_compound_id],
        namespace="kegg.compound",
    )


def _chebi_matches(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _chebi_selection(inputs)
    return {
        "matches": _collect(
            selection.matches(),
            "input_id",
            "input_namespace",
            "chebi_id",
        )
    }


def _chebi_bundle(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _chebi_selection(inputs)
    return {
        "matches": _collect(
            selection.matches(),
            "input_id",
            "input_namespace",
            "chebi_id",
        ),
        "compounds": _collect(
            selection.compounds(),
            "input_id",
            "input_namespace",
            "chebi_id",
        ),
        "unmatched": _collect(
            selection.unmatched_ids(),
            "input_id",
            "input_namespace",
            "reason",
        ),
    }


def _kegg_selection(inputs: BenchmarkInputs):
    return KEGGDatabase.from_duckdb(inputs.kegg).select_ids(
        [inputs.kegg_compound_id],
        namespace="kegg_compound",
    )


def _kegg_matches(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _kegg_selection(inputs)
    return {"matches": _collect(selection.matches(), "input_id")}


def _kegg_reactions_narrow(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _kegg_selection(inputs)
    return {
        "reactions": _collect(
            selection.reactions(),
            "input_id",
            "reaction_id",
            "is_reversible",
        )
    }


def _kegg_reactions_wide(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _kegg_selection(inputs)
    return {"reactions": selection.reactions().collect(engine="streaming")}


def _kegg_bundle(inputs: BenchmarkInputs) -> Mapping[str, pl.DataFrame]:
    selection = _kegg_selection(inputs)
    return {
        "matches": _collect(selection.matches(), "input_id"),
        "reactions": _collect(
            selection.reactions(),
            "input_id",
            "reaction_id",
            "is_reversible",
        ),
        "participants": _collect(
            selection.participants(),
            "input_id",
            "reaction_id",
            "side",
            "participant_id",
            "coefficient_numeric",
        ),
        "enzymes": _collect(
            selection.enzymes(),
            "input_id",
            "reaction_id",
            "ec_number",
        ),
        "kos": _collect(
            selection.kos(),
            "input_id",
            "reaction_id",
            "ko_id",
        ),
        "modules": _collect(
            selection.modules(),
            "input_id",
            "reaction_id",
            "module_id",
        ),
        "pathways": _collect(
            selection.pathway_memberships(),
            "input_id",
            "reaction_id",
            "pathway_id",
        ),
        "unmatched": _collect(selection.unmatched_ids(), "input_id", "reason"),
    }


def _case_operations(inputs: BenchmarkInputs) -> Mapping[str, CaseOperation]:
    functions: Mapping[
        str,
        Callable[[BenchmarkInputs], Mapping[str, pl.DataFrame]],
    ] = {
        "idmapping_narrow": _idmapping_narrow,
        "idmapping_wide": _idmapping_wide,
        "rhea_matches": _rhea_matches,
        "rhea_bundle": _rhea_bundle,
        "uniprot_proteins": _uniprot_proteins,
        "uniprot_bundle": _uniprot_bundle,
        "chebi_matches": _chebi_matches,
        "chebi_bundle": _chebi_bundle,
        "kegg_matches": _kegg_matches,
        "kegg_reactions_narrow": _kegg_reactions_narrow,
        "kegg_reactions_wide": _kegg_reactions_wide,
        "kegg_bundle": _kegg_bundle,
    }
    return {name: lambda fn=fn: fn(inputs) for name, fn in functions.items()}


def _json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        default=str,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()


def _sha256(value: object) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _frame_contract(frame: pl.DataFrame) -> dict[str, object]:
    schema = [(name, str(dtype)) for name, dtype in frame.schema.items()]
    records = frame.to_dicts()
    return {
        "rows": frame.height,
        "columns": frame.columns,
        "schema": schema,
        "schema_sha256": _sha256(schema),
        "semantic_sha256": _sha256(records),
        "estimated_size_bytes": frame.estimated_size(),
    }


def _output_contract(frames: Mapping[str, pl.DataFrame]) -> dict[str, object]:
    relations = {name: _frame_contract(frame) for name, frame in sorted(frames.items())}
    return {
        "relations": relations,
        "contract_sha256": _sha256(relations),
    }


def _measure(operation: CaseOperation) -> dict[str, object]:
    gc.collect()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    wall_started = time.perf_counter()
    frames = operation()
    wall_seconds = time.perf_counter() - wall_started
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    contract = _output_contract(frames)
    return {
        "wall_s": wall_seconds,
        "user_s": usage_after.ru_utime - usage_before.ru_utime,
        "system_s": usage_after.ru_stime - usage_before.ru_stime,
        "cpu_s": (
            usage_after.ru_utime
            - usage_before.ru_utime
            + usage_after.ru_stime
            - usage_before.ru_stime
        ),
        "peak_rss_kib_process_global": usage_after.ru_maxrss,
        **contract,
    }


def _metric_summary(
    samples: Sequence[Mapping[str, object]],
    metric: str,
) -> dict[str, float]:
    values: list[float] = []
    for sample in samples:
        raw_value = sample[metric]
        if not isinstance(raw_value, int | float):
            raise TypeError(f"benchmark metric {metric!r} is not numeric")
        values.append(float(raw_value))
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def _aggregate(samples: Sequence[Mapping[str, object]]) -> dict[str, object]:
    hashes = {str(sample["contract_sha256"]) for sample in samples}
    if len(hashes) != 1:
        raise RuntimeError(f"case output changed across runs: {sorted(hashes)}")
    return {
        "wall_s": _metric_summary(samples, "wall_s"),
        "cpu_s": _metric_summary(samples, "cpu_s"),
        "user_s": _metric_summary(samples, "user_s"),
        "system_s": _metric_summary(samples, "system_s"),
        "peak_rss_kib_process_global": _metric_summary(
            samples,
            "peak_rss_kib_process_global",
        ),
        "contract_sha256": hashes.pop(),
        "relations": samples[0]["relations"],
    }


def _publication_identity(path: Path) -> dict[str, object]:
    inspection = inspect_publication(path)
    stat = path.stat()
    contract = {
        "metadata": [asdict(record) for record in inspection.metadata],
        "tables": [asdict(record) for record in inspection.tables],
        "source_files": [asdict(record) for record in inspection.source_files],
        "column_mappings": [asdict(record) for record in inspection.column_mappings],
        "validation_issues": [
            asdict(record) for record in inspection.validation_issues
        ],
    }
    return {
        "path": str(path.resolve()),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "resource_name": inspection.resource_name,
        "resource_schema_version": inspection.resource_schema_version,
        "source_schema_profile": inspection.source_schema_profile,
        "release_version": inspection.release_version,
        "package_version": inspection.package_version,
        "generated_at": inspection.generated_at,
        "scope": inspection.scope,
        "validation_status": inspection.validation_status,
        "validation_issue_count": inspection.validation_issue_count,
        "table_count": len(inspection.tables),
        "recorded_rows": sum(table.row_count for table in inspection.tables),
        "publication_contract_sha256": _sha256(contract),
    }


def _duckdb_threads() -> int:
    with duckdb.connect(":memory:") as connection:
        value = connection.execute("SELECT current_setting('threads')").fetchone()
    if value is None:
        raise RuntimeError("DuckDB did not report its thread setting")
    return int(value[0])


def _package_version() -> str:
    try:
        return importlib.metadata.version("bioextract")
    except importlib.metadata.PackageNotFoundError:
        return "working-tree"


def _environment(label: str) -> dict[str, object]:
    return {
        "label": label,
        "benchmark_script_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "bioextract": _package_version(),
        "polars": pl.__version__,
        "duckdb": duckdb.__version__,
        "polars_thread_pool_size": pl.thread_pool_size(),
        "duckdb_threads": _duckdb_threads(),
        "cpu_count": os.cpu_count(),
        "thread_environment": {
            name: os.environ[name] for name in _THREAD_ENVIRONMENT if name in os.environ
        },
    }


def _run(
    *,
    inputs: BenchmarkInputs,
    cases: Sequence[str],
    runs: int,
    warmups: int,
    label: str,
) -> dict[str, object]:
    operations = _case_operations(inputs)
    for warmup in range(warmups):
        for case in cases:
            print(
                f"warmup {warmup + 1}/{warmups}: {case}",
                file=sys.stderr,
                flush=True,
            )
            operations[case]()

    samples_by_case: dict[str, list[dict[str, object]]] = {case: [] for case in cases}
    for run_index in range(runs):
        ordered_cases = list(cases)
        if run_index % 2:
            ordered_cases.reverse()
        for order, case in enumerate(ordered_cases, 1):
            print(
                f"run {run_index + 1}/{runs}, {order}/{len(cases)}: {case}",
                file=sys.stderr,
                flush=True,
            )
            sample = _measure(operations[case])
            sample.update({"run": run_index + 1, "order": order})
            samples_by_case[case].append(sample)

    paths = {
        "idmapping": inputs.idmapping,
        "uniprot_kb": inputs.uniprot_kb,
        "rhea": inputs.rhea,
        "chebi": inputs.chebi,
        "kegg": inputs.kegg,
    }
    return {
        "schema_version": "bioextract-lazy-relation-benchmark-v1",
        "environment": _environment(label),
        "inputs": {
            "taxon_id": inputs.taxon_id,
            "protein_ids": inputs.protein_ids,
            "kegg_compound_id": inputs.kegg_compound_id,
            "publications": {
                name: _publication_identity(path) for name, path in paths.items()
            },
        },
        "protocol": {
            "runs": runs,
            "warmups": warmups,
            "order": "forward on odd runs and reverse on even runs",
            "cases": list(cases),
        },
        "summary": {case: _aggregate(samples_by_case[case]) for case in cases},
        "samples": samples_by_case,
        "notes": [
            "relations use only public bioextract APIs",
            "each timed case reconstructs its database handle and selection",
            "peak RSS is process-global and is not isolated per case",
            "publication inspection is outside timed cases",
        ],
    }


def main() -> None:
    args = _parser().parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be at least 1")
    if args.warmups < 0:
        raise SystemExit("--warmups must be non-negative")
    protein_ids = tuple(dict.fromkeys(str(value) for value in args.protein_id))
    if not protein_ids:
        raise SystemExit("at least one non-empty --protein-id is required")
    cases = tuple(dict.fromkeys(args.case or _CASE_NAMES))
    inputs = BenchmarkInputs(
        idmapping=args.idmapping,
        uniprot_kb=args.uniprot_kb,
        rhea=args.rhea,
        chebi=args.chebi,
        kegg=args.kegg,
        taxon_id=args.taxon_id,
        protein_ids=protein_ids,
        kegg_compound_id=args.kegg_compound_id,
    )
    report = _run(
        inputs=inputs,
        cases=cases,
        runs=args.runs,
        warmups=args.warmups,
        label=args.label,
    )
    rendered = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(f"{rendered}\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
