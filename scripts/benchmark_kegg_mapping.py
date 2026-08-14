"""Benchmark one KEGG mapping publication with bounded observable evidence."""

from __future__ import annotations

import argparse
import cProfile
import hashlib
import json
import os
import resource
import threading
import time
from collections.abc import Callable, Generator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import duckdb

from bioextract import KEGGDatabase
from bioextract.kegg.mapping import _native as native_mapping
from bioextract.kegg.mapping import publication as mapping_publication


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--organism",
        action="append",
        default=[],
        help="Organism code to publish; repeat to benchmark a bounded subset.",
    )
    parser.add_argument(
        "--organism-count",
        type=int,
        help="Benchmark the first N lexically sorted organism directories.",
    )
    parser.add_argument(
        "--organism-file",
        type=Path,
        help="Read the exact organism scope from a newline-delimited file.",
    )
    parser.add_argument(
        "--expect-organism-scope-sha256",
        help="Fail before building unless the normalized scope has this SHA-256.",
    )
    parser.add_argument(
        "--profile-output",
        type=Path,
        help="Write cProfile data for write_duckdb(); profiling changes timing.",
    )
    parser.add_argument(
        "--phase-report",
        action="store_true",
        help="Record internal writer phase timings; adds wrapper overhead.",
    )
    parser.add_argument("--release-version")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    organism_scope = list(args.organism)
    if args.organism_file is not None:
        if organism_scope or args.organism_count is not None:
            parser.error(
                "--organism-file, --organism, and --organism-count are "
                "mutually exclusive"
            )
        organism_scope = [
            line.strip()
            for line in args.organism_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not organism_scope:
            parser.error("--organism-file must contain at least one organism code")
    if args.organism_count is not None:
        if args.organism_count <= 0:
            parser.error("--organism-count must be positive")
        if organism_scope:
            parser.error("--organism and --organism-count are mutually exclusive")
        organism_scope = sorted(
            entry.name
            for entry in os.scandir(args.source)
            if entry.is_dir(follow_symlinks=False)
            and entry.name not in {"ko", "organism"}
            and 3 <= len(entry.name) <= 4
            and entry.name.islower()
            and entry.name.isalpha()
        )[: args.organism_count]
    organism_scope = sorted(set(organism_scope))
    scope_sha256 = (
        hashlib.sha256("\n".join(organism_scope).encode()).hexdigest()
        if organism_scope
        else None
    )
    if (
        args.expect_organism_scope_sha256 is not None
        and scope_sha256 != args.expect_organism_scope_sha256
    ):
        parser.error(
            "normalized organism scope SHA-256 mismatch: "
            f"expected {args.expect_organism_scope_sha256}, observed {scope_sha256}"
        )
    database = KEGGDatabase.from_mapping_directory(
        args.source,
        release_version=args.release_version,
    )
    if organism_scope:
        database = database.with_organisms(organism_scope)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    samples: list[dict[str, float | int]] = []
    memory_limit, max_temp_directory_size = native_mapping._resource_limits(  # pyright: ignore[reportPrivateUsage]
        args.output.parent
    )
    duckdb_threads = min(
        native_mapping.pl.thread_pool_size(),
        native_mapping._THREAD_CAP,  # pyright: ignore[reportPrivateUsage]
    )
    stop_sampling = threading.Event()
    monitor = threading.Thread(
        target=_sample_process,
        args=(stop_sampling, samples, args.output.parent, args.output.name),
        daemon=True,
    )
    monitor.start()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    started = time.perf_counter()
    profile_output: Path | None = args.profile_output
    profiler = cProfile.Profile() if profile_output is not None else None
    with _phase_instrumentation(args.phase_report) as phase_seconds:
        try:
            if profiler is not None:
                profiler.enable()
            try:
                result = database.write_duckdb(
                    args.output,
                    if_exists="replace" if args.replace else "fail",
                )
            finally:
                if profiler is not None:
                    profiler.disable()
                    if profile_output is None:
                        raise AssertionError("profile output disappeared")
                    profile_output.parent.mkdir(parents=True, exist_ok=True)
                    profiler.dump_stats(profile_output)
        finally:
            stop_sampling.set()
            monitor.join()
    elapsed = time.perf_counter() - started
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    with KEGGDatabase.from_duckdb(result.path).connect() as connection:
        source_summary = connection.execute(
            "SELECT count(*), count(bytes), sum(bytes) FROM _bioextract.source_file"
        ).fetchone()
        metadata = dict(
            connection.execute("SELECT key, value FROM _bioextract.metadata").fetchall()
        )
    if source_summary is None:
        raise RuntimeError("KEGG benchmark could not read source inventory")
    report: dict[str, Any] = {
        "source": str(args.source),
        "output": str(result.path),
        "organism_scope": (
            sorted(set(organism_scope)) if len(organism_scope) <= 20 else None
        ),
        "organism_scope_count": len(set(organism_scope)) or None,
        "organism_scope_sha256": scope_sha256,
        "elapsed_seconds": elapsed,
        "cpu_seconds": (
            usage_after.ru_utime
            + usage_after.ru_stime
            - usage_before.ru_utime
            - usage_before.ru_stime
        ),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "sampled_peak_rss_bytes": max(
            (int(sample["rss_bytes"]) for sample in samples), default=0
        ),
        "sampled_peak_open_file_descriptors": max(
            (int(sample["open_file_descriptors"]) for sample in samples), default=0
        ),
        "sampled_peak_temp_bytes": max(
            (int(sample["temp_bytes"]) for sample in samples), default=0
        ),
        "sampled_peak_stage_bytes": max(
            (int(sample["stage_bytes"]) for sample in samples), default=0
        ),
        "source_file_count": int(source_summary[0]),
        "source_bytes": (
            None
            if int(source_summary[1]) != int(source_summary[0])
            else int(source_summary[2] or 0)
        ),
        "source_bytes_known": int(source_summary[1]),
        "output_bytes": result.path.stat().st_size,
        "row_counts": dict(result.row_counts),
        "validation_issue_count": result.validation_issue_count,
        "profiling_enabled": profiler is not None,
        "phase_seconds": phase_seconds,
        "output_sha256": _sha256_file(result.path),
        "capabilities": {
            key.removeprefix("bioextract.capability."): value == "true"
            for key, value in metadata.items()
            if key.startswith("bioextract.capability.")
        },
        "duckdb": {
            "version": duckdb.__version__,
            "threads": duckdb_threads,
            "memory_limit": memory_limit,
            "max_temp_directory_size": max_temp_directory_size,
            "temp_directory": "per-publication-staged-temporary-directory",
        },
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")


def _sample_process(
    stop: threading.Event,
    samples: list[dict[str, float | int]],
    output_parent: Path,
    output_name: str,
) -> None:
    started = time.monotonic()
    while not stop.wait(1.0):
        rss_bytes = 0
        try:
            for line in Path("/proc/self/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    rss_bytes = int(line.split()[1]) * 1024
                    break
            descriptor_count = sum(1 for _ in Path("/proc/self/fd").iterdir())
            temp_bytes = _matching_directory_bytes(
                output_parent, "bioextract-kegg-duckdb-"
            )
            stage_bytes = _matching_file_bytes(output_parent, f".{output_name}.")
        except OSError:
            descriptor_count = 0
            temp_bytes = 0
            stage_bytes = 0
        samples.append(
            {
                "elapsed_seconds": time.monotonic() - started,
                "rss_bytes": rss_bytes,
                "open_file_descriptors": descriptor_count,
                "temp_bytes": temp_bytes,
                "stage_bytes": stage_bytes,
            }
        )


def _matching_directory_bytes(parent: Path, prefix: str) -> int:
    total = 0
    for candidate in parent.iterdir():
        if not candidate.is_dir() or not candidate.name.startswith(prefix):
            continue
        for path in candidate.rglob("*"):
            if path.is_file():
                total += path.stat().st_size
    return total


def _matching_file_bytes(parent: Path, prefix: str) -> int:
    return sum(
        path.stat().st_size
        for path in parent.iterdir()
        if path.is_file() and path.name.startswith(prefix)
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def _phase_instrumentation(
    enabled: bool,
) -> Generator[dict[str, float]]:
    """Temporarily time private KEGG writer phases for benchmark evidence."""
    phases: dict[str, float] = {}
    if not enabled:
        yield phases
        return

    originals: list[tuple[Any, str, Callable[..., Any]]] = []

    def wrap(module: Any, name: str, phase: str) -> None:
        original = getattr(module, name)

        def timed(*args: Any, **kwargs: Any) -> Any:
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                phases[phase] = phases.get(phase, 0.0) + (time.perf_counter() - started)

        originals.append((module, name, original))
        setattr(module, name, timed)

    wrap(native_mapping, "_scan_role", "scan_role")
    wrap(native_mapping, "_validate_role_content", "role_validation")
    wrap(native_mapping, "_build_gene_table", "gene_table")
    wrap(native_mapping, "_build_ko_table", "ko_table")
    wrap(mapping_publication, "validate_mapping_publication", "mapping_validation")
    try:
        yield phases
    finally:
        for module, name, original in reversed(originals):
            setattr(module, name, original)


if __name__ == "__main__":
    main()
