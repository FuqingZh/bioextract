"""Evaluate a native multi-file KEGG role scan without changing the publisher."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import polars as pl
from polars._typing import SchemaDict

from bioextract.kegg.mapping import parse as mapping_parse
from bioextract.kegg.mapping.source import validate_organism_code

ROLE_COLUMNS: dict[str, tuple[str, ...]] = {
    "gene_list": ("kegg_gene_id", "gene_type", "genomic_position", "display"),
    "uniprot_conversion": ("xref", "kegg_gene_id"),
    "ncbi_gene_conversion": ("xref", "kegg_gene_id"),
    "gene_ko": ("kegg_gene_id", "ko_id"),
    "gene_pathway": ("kegg_gene_id", "pathway_id"),
}
ROLE_FILENAMES: dict[str, str] = {
    "gene_list": "gene_list.tsv",
    "uniprot_conversion": "conv_uniprot.tsv",
    "ncbi_gene_conversion": "conv_ncbi_geneid.tsv",
    "gene_ko": "gene_ko.tsv",
    "gene_pathway": "gene_pathway.tsv",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--engine",
        choices=("python", "polars"),
        required=True,
        help="Read roles with the current parser or the native Polars prototype.",
    )
    parser.add_argument(
        "--organism-count",
        type=int,
        help="Use the first N lexically sorted organism directories.",
    )
    parser.add_argument(
        "--organism-file",
        type=Path,
        help="Read an exact newline-delimited organism scope.",
    )
    parser.add_argument("--expect-scope-sha256")
    parser.add_argument("--report", type=Path)
    parser.add_argument(
        "--skip-contract-probe",
        action="store_true",
        help="Skip the temporary malformed/empty-file Polars probe.",
    )
    args = parser.parse_args()

    codes = _read_scope(
        args.source,
        organism_count=args.organism_count,
        organism_file=args.organism_file,
    )
    scope_sha256 = hashlib.sha256("\n".join(codes).encode()).hexdigest()
    if (
        args.expect_scope_sha256 is not None
        and scope_sha256 != args.expect_scope_sha256
    ):
        parser.error(
            "normalized organism scope SHA-256 mismatch: "
            f"expected {args.expect_scope_sha256}, observed {scope_sha256}"
        )
    role_paths = _role_paths(args.source, codes)
    started = time.perf_counter()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    if args.engine == "python":
        row_counts = _read_python(role_paths)
        native_details: dict[str, Any] = {}
    else:
        row_counts, native_details = _read_polars(role_paths)
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    report: dict[str, Any] = {
        "engine": args.engine,
        "source": str(args.source),
        "scope_count": len(codes),
        "scope_sha256": scope_sha256,
        "source_file_count": sum(len(paths) for paths in role_paths.values()),
        "enumerated_source_bytes": sum(
            path.stat().st_size for paths in role_paths.values() for path in paths
        ),
        "elapsed_seconds": time.perf_counter() - started,
        "cpu_seconds": (
            usage_after.ru_utime
            + usage_after.ru_stime
            - usage_before.ru_utime
            - usage_before.ru_stime
        ),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "row_counts": row_counts,
        "role_file_counts": {role: len(paths) for role, paths in role_paths.items()},
        "native_details": native_details,
        "environment": {
            "python": platform.python_version(),
            "polars": pl.__version__,
        },
        "contract_probe": (None if args.skip_contract_probe else _contract_probe()),
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")


def _read_scope(
    source: Path,
    *,
    organism_count: int | None,
    organism_file: Path | None,
) -> tuple[str, ...]:
    if organism_count is not None and organism_file is not None:
        raise ValueError("--organism-count and --organism-file are mutually exclusive")
    if organism_file is not None:
        codes = tuple(
            sorted(
                {
                    validate_organism_code(line.strip())
                    for line in organism_file.read_text(encoding="utf-8").splitlines()
                    if line.strip() and not line.lstrip().startswith("#")
                }
            )
        )
    else:
        if organism_count is None:
            raise ValueError("one of --organism-count or --organism-file is required")
        if organism_count <= 0:
            raise ValueError("--organism-count must be positive")
        with os.scandir(source) as entries:
            codes = tuple(
                sorted(
                    entry.name
                    for entry in entries
                    if entry.is_dir(follow_symlinks=False)
                    and entry.name not in {"ko", "organism"}
                    and 3 <= len(entry.name) <= 4
                    and entry.name.islower()
                    and entry.name.isalpha()
                )[:organism_count]
            )
    if not codes:
        raise ValueError("organism scope must contain at least one code")
    return codes


def _role_paths(source: Path, codes: tuple[str, ...]) -> dict[str, list[Path]]:
    paths: dict[str, list[Path]] = {role: [] for role in ROLE_FILENAMES}
    for code in codes:
        for role, filename in ROLE_FILENAMES.items():
            path = source / code / filename
            if not path.is_file():
                raise FileNotFoundError(path)
            paths[role].append(path)
    return paths


def _read_python(role_paths: dict[str, list[Path]]) -> dict[str, int]:
    row_counts: dict[str, int] = {}
    for role, paths in role_paths.items():
        rows = 0
        for path in paths:
            parsed = mapping_parse._read_rows(  # pyright: ignore[reportPrivateUsage]
                path,
                logical_name=f"p3/{role}/{path.name}",
                columns=len(ROLE_COLUMNS[role]),
            )
            rows += len(parsed.rows)
        row_counts[role] = rows
    return row_counts


def _read_polars(
    role_paths: dict[str, list[Path]],
) -> tuple[dict[str, int], dict[str, Any]]:
    row_counts: dict[str, int] = {}
    path_counts: dict[str, int] = {}
    for role, paths in role_paths.items():
        columns = ROLE_COLUMNS[role]
        schema: SchemaDict = dict.fromkeys(columns, pl.String)
        lf = pl.scan_csv(
            paths,
            has_header=False,
            separator="\t",
            schema=schema,
            include_file_paths="_source_path",
            raise_if_empty=False,
            truncate_ragged_lines=False,
            ignore_errors=False,
            low_memory=True,
            rechunk=False,
            glob=False,
        ).with_columns(pl.col(name).str.strip_chars() for name in columns)
        frame = lf.unique(subset=["_source_path", *columns]).collect()
        path_to_organism = {str(path): path.parent.name for path in paths}
        mapped = frame.with_columns(
            pl.col("_source_path").replace(path_to_organism).alias("_organism_code")
        )
        if mapped.get_column("_organism_code").null_count():
            raise AssertionError(f"unmapped source path in role {role!r}")
        row_counts[role] = frame.height
        path_counts[role] = frame.get_column("_source_path").n_unique()
    return row_counts, {
        "per_file_unique": True,
        "same_read_source_hash": False,
        "path_counts": path_counts,
    }


def _contract_probe() -> dict[str, Any]:
    columns = ("left", "right")
    schema: SchemaDict = dict.fromkeys(columns, pl.String)

    def scan(path: Path) -> pl.DataFrame:
        return pl.scan_csv(
            [path],
            has_header=False,
            separator="\t",
            schema=schema,
            include_file_paths="_source_path",
            raise_if_empty=False,
            truncate_ragged_lines=False,
            ignore_errors=False,
            glob=False,
        ).collect()

    with TemporaryDirectory() as directory:
        root = Path(directory)
        empty = root / "empty.tsv"
        good = root / "good.tsv"
        missing = root / "missing.tsv"
        extra = root / "extra.tsv"
        invalid_utf8 = root / "invalid-utf8.tsv"
        empty.write_bytes(b"")
        good.write_text("a\tb\n", encoding="utf-8")
        missing.write_text("a\n", encoding="utf-8")
        extra.write_text("a\tb\tc\n", encoding="utf-8")
        invalid_utf8.write_bytes(b"a\t\xff\n")
        empty_frame = scan(empty)
        good_empty = pl.scan_csv(
            [good, empty],
            has_header=False,
            separator="\t",
            schema=schema,
            include_file_paths="_source_path",
            raise_if_empty=False,
            truncate_ragged_lines=False,
            ignore_errors=False,
            glob=False,
        ).collect()
        errors = {
            path.name: _capture_scan_error(scan, path)
            for path in (missing, extra, invalid_utf8)
        }
        return {
            "zero_byte_collects": empty_frame.height == 0,
            "zero_byte_appears_in_path_column": str(empty)
            in set(good_empty.get_column("_source_path").to_list()),
            "good_plus_zero_byte_rows": good_empty.height,
            "errors": errors,
            "malformed_error_includes_path": any(
                str(path) in error["message"]
                for path, error in (
                    (missing, errors[missing.name]),
                    (extra, errors[extra.name]),
                    (invalid_utf8, errors[invalid_utf8.name]),
                )
            ),
            "malformed_error_includes_line": any(
                "line" in error["message"].lower() for error in errors.values()
            ),
        }


def _capture_scan_error(scan: Any, path: Path) -> dict[str, str]:
    try:
        scan(path)
    except Exception as error:
        return {"type": type(error).__name__, "message": str(error)}
    raise AssertionError(f"expected native scan failure for {path}")


if __name__ == "__main__":
    main()
