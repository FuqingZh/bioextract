from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
import tomllib
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePosixPath
from typing import Any, cast

import tomli_w

from bioextract.publication import (
    PublicationDescriptor,
    inspect_duckdb_publication,
)

CATALOG_SCHEMA_VERSION = "bioextract-publication-catalog-v1"
BIOFETCH_MANIFEST_SCHEMA_VERSION = "biofetch-manifest-v1"
PUBLICATION_RELATIVE_PATH = PurePosixPath("tidy/data.duckdb")


def build_publication_catalog(
    manifest_path: os.PathLike[str] | str,
    output_path: os.PathLike[str] | str,
) -> dict[str, Any]:
    """Build one deterministic catalog of validated DuckDB publications."""
    source_manifest_path = Path(manifest_path).resolve(strict=True)
    if not source_manifest_path.is_file():
        raise IsADirectoryError(source_manifest_path)
    destination = Path(output_path).resolve()
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)

    manifest_bytes = source_manifest_path.read_bytes()
    manifest = tomllib.loads(manifest_bytes.decode("utf-8"))
    _require_string(
        manifest,
        "schema_version",
        expected=BIOFETCH_MANIFEST_SCHEMA_VERSION,
    )
    resource_root_value = _require_string(manifest, "resource_root")
    resource_root = (source_manifest_path.parent / resource_root_value).resolve(
        strict=True
    )
    if not resource_root.is_dir():
        raise NotADirectoryError(resource_root)

    snapshots = _require_list(manifest, "snapshots")
    publications: list[dict[str, Any]] = []
    identities: set[tuple[str, str, str]] = set()
    for index, snapshot_value in enumerate(snapshots):
        snapshot = _require_mapping(snapshot_value, f"snapshots[{index}]")
        database = _require_string(snapshot, "database")
        asset = _require_string(snapshot, "asset")
        version = _require_string(snapshot, "version")
        identity = (database, asset, version)
        if identity in identities:
            raise ValueError(f"Duplicate snapshot identity: {identity}")
        identities.add(identity)

        snapshot_relative = _require_relative_path(
            _require_string(snapshot, "path"),
            context=f"snapshots[{index}].path",
        )
        snapshot_root = (resource_root / snapshot_relative).resolve(strict=True)
        _require_within(snapshot_root, resource_root, context="snapshot")
        if not snapshot_root.is_dir():
            raise NotADirectoryError(snapshot_root)

        lock = _require_mapping(
            snapshot.get("manifest"),
            f"snapshots[{index}].manifest",
        )
        lock_relative = _require_relative_path(
            _require_string(lock, "path"),
            context=f"snapshots[{index}].manifest.path",
        )
        lock_path = (snapshot_root / lock_relative).resolve(strict=True)
        _require_within(lock_path, snapshot_root, context="snapshot manifest")
        _validate_declared_file(
            lock_path,
            declared_sha256=_require_string(lock, "sha256"),
            declared_bytes=_require_integer(lock, "bytes"),
            context="snapshot manifest",
        )

        publication_candidate = snapshot_root / PUBLICATION_RELATIVE_PATH
        if not publication_candidate.exists():
            continue
        publication_path = publication_candidate.resolve(strict=True)
        _require_within(
            publication_path,
            snapshot_root,
            context="DuckDB publication",
        )
        descriptor = inspect_duckdb_publication(publication_path)
        publications.append(
            _publication_record(
                descriptor,
                database=database,
                asset=asset,
                version=version,
                snapshot_path=snapshot_relative.as_posix(),
                lock_path=lock_relative.as_posix(),
                lock_sha256=_require_string(lock, "sha256"),
                lock_bytes=_require_integer(lock, "bytes"),
            )
        )

    publications.sort(
        key=lambda item: (
            str(item["database"]),
            str(item["asset"]),
            str(item["version"]),
        )
    )
    relative_root = Path(os.path.relpath(resource_root, destination.parent)).as_posix()
    relative_manifest = Path(
        os.path.relpath(source_manifest_path, destination.parent)
    ).as_posix()
    model: dict[str, Any] = {
        "schema_version": CATALOG_SCHEMA_VERSION,
        "resource_root": relative_root,
        "source_manifest": {
            "path": relative_manifest,
            "schema_version": BIOFETCH_MANIFEST_SCHEMA_VERSION,
            "sha256": _sha256_bytes(manifest_bytes),
            "bytes": len(manifest_bytes),
        },
        "summary": {
            "snapshot_count": len(snapshots),
            "publication_count": len(publications),
            "total_bytes": sum(int(item["bytes"]) for item in publications),
        },
        "publications": publications,
    }
    _write_toml_atomically(destination, model)
    return model


def _publication_record(
    descriptor: PublicationDescriptor,
    *,
    database: str,
    asset: str,
    version: str,
    snapshot_path: str,
    lock_path: str,
    lock_sha256: str,
    lock_bytes: int,
) -> dict[str, Any]:
    publication: dict[str, Any] = {
        "database": database,
        "asset": asset,
        "version": version,
        "snapshot_path": snapshot_path,
        "path": PUBLICATION_RELATIVE_PATH.as_posix(),
        "sha256": _sha256_file(descriptor.path),
        "bytes": descriptor.path.stat().st_size,
        "source_manifest_path": lock_path,
        "source_manifest_sha256": lock_sha256,
        "source_manifest_bytes": lock_bytes,
        "metadata_schema_version": descriptor.metadata_schema_version,
        "resource_name": descriptor.resource_name,
        "resource_schema_version": descriptor.resource_schema_version,
        "source_schema_profile": descriptor.source_schema_profile,
        "package_version": descriptor.package_version,
        "generated_at": descriptor.generated_at,
        "validation_status": descriptor.validation_status,
        "validation_issue_count": descriptor.validation_issue_count,
        "source_count": len(descriptor.sources),
        "tables": [
            {
                "name": table.name,
                "role": table.role,
                "row_count": table.row_count,
            }
            for table in descriptor.tables
        ],
    }
    optional_values = {
        "source_schema_version": descriptor.source_schema_version,
        "release_version": descriptor.release_version,
        "release_version_source": descriptor.release_version_source,
        "scope": descriptor.scope,
    }
    publication.update(
        {key: value for key, value in optional_values.items() if value is not None}
    )
    return publication


def _validate_declared_file(
    path: Path,
    *,
    declared_sha256: str,
    declared_bytes: int,
    context: str,
) -> None:
    if not path.is_file():
        raise IsADirectoryError(path)
    observed_bytes = path.stat().st_size
    if observed_bytes != declared_bytes:
        raise ValueError(
            f"{context} byte count does not match: "
            f"declared={declared_bytes}, actual={observed_bytes}"
        )
    observed_sha256 = _sha256_file(path)
    if observed_sha256 != declared_sha256:
        raise ValueError(
            f"{context} SHA-256 does not match: "
            f"declared={declared_sha256}, actual={observed_sha256}"
        )


def _write_toml_atomically(destination: Path, model: Mapping[str, Any]) -> None:
    payload = tomli_w.dumps(dict(model), multiline_strings=False).encode("utf-8")
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        delete=False,
    ) as handle:
        stage = Path(handle.name)
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(stage, destination)
    finally:
        stage.unlink(missing_ok=True)


def _require_mapping(value: object, context: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a TOML table")
    untyped_mapping = cast(dict[object, object], value)
    if not all(isinstance(key, str) for key in untyped_mapping):
        raise ValueError(f"{context} must be a TOML table")
    return cast(dict[str, Any], untyped_mapping)


def _require_list(mapping: Mapping[str, Any], key: str) -> list[object]:
    value = mapping.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be an array")
    return cast(list[object], value)


def _require_string(
    mapping: Mapping[str, Any],
    key: str,
    *,
    expected: str | None = None,
) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{key} must be a non-empty string")
    if expected is not None and value != expected:
        raise ValueError(f"{key} must be {expected!r}, got {value!r}")
    return value


def _require_integer(mapping: Mapping[str, Any], key: str) -> int:
    value = mapping.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{key} must be a non-negative integer")
    return value


def _require_relative_path(value: str, *, context: str) -> PurePosixPath:
    if "\\" in value:
        raise ValueError(f"{context} must use forward slashes: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or value != path.as_posix() or ".." in path.parts:
        raise ValueError(f"{context} must be a canonical relative path: {value!r}")
    return path


def _require_within(path: Path, root: Path, *, context: str) -> None:
    if not path.is_relative_to(root):
        raise ValueError(f"{context} resolves outside its owning root: {path}")


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a deterministic TOML catalog from biofetch snapshots that "
            "contain validated tidy/data.duckdb publications."
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="biofetch aggregate manifest.toml",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="destination catalog.toml; its parent directory must exist",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    model = build_publication_catalog(args.manifest, args.output)
    summary = _require_mapping(model["summary"], "summary")
    print(
        "Built publication catalog: "
        f"snapshots={summary['snapshot_count']} "
        f"publications={summary['publication_count']} "
        f"output={args.output.resolve()}"
    )


if __name__ == "__main__":
    main()
