from __future__ import annotations

import hashlib
import runpy
import tomllib
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, cast

import duckdb
import polars as pl
import pytest
import tomli_w

from bioextract._tidy import TidyAsset, TidyDataset, TidySource

_SCRIPT_GLOBALS = runpy.run_path(
    str(Path(__file__).parents[1] / "scripts/build_publication_catalog.py")
)
build_publication_catalog = cast(
    Callable[[Path, Path], dict[str, Any]],
    _SCRIPT_GLOBALS["build_publication_catalog"],
)
main = cast(Callable[[Sequence[str] | None], None], _SCRIPT_GLOBALS["main"])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot(
    resource_root: Path,
    *,
    database: str,
    asset: str,
    version: str,
    publication: bool,
) -> dict[str, object]:
    relative = Path(database) / asset / version
    snapshot = resource_root / relative
    raw = snapshot / "raw"
    raw.mkdir(parents=True)
    source = raw / "source.tsv"
    source.write_text("id\nT1\n", encoding="utf-8")
    lock = snapshot / "manifest.lock"
    lock.write_text(
        f"database = {database!r}\nasset = {asset!r}\nversion_token = {version!r}\n",
        encoding="utf-8",
    )
    if publication:
        tidy = snapshot / "tidy"
        tidy.mkdir()
        TidyDataset(
            frames={"mapping": pl.DataFrame({"id": ["T1"]}).lazy()},
            source=TidySource("source", source, "text/tab-separated-values"),
            resource_schema_version=f"{database}-v1",
            source_schema_profile=f"{database}-{asset}-v1",
            build_id_prefix=database,
            assets=(TidyAsset("mapping.parquet", "canonical", "mapping"),),
            resource_name=database,
            release_version=version,
        ).write_duckdb(tidy / "data.duckdb")
    return {
        "database": database,
        "asset": asset,
        "version": version,
        "path": relative.as_posix(),
        "manifest": {
            "path": "manifest.lock",
            "sha256": _sha256(lock),
            "bytes": lock.stat().st_size,
        },
    }


def _write_manifest(
    tmp_path: Path,
    *,
    snapshots: list[dict[str, object]],
) -> tuple[Path, Path]:
    meta = tmp_path / "meta"
    meta.mkdir()
    manifest = meta / "manifest.toml"
    manifest.write_bytes(
        tomli_w.dumps(
            {
                "schema_version": "biofetch-manifest-v1",
                "resource_root": "../resources",
                "snapshots": snapshots,
            }
        ).encode("utf-8")
    )
    return manifest, meta / "catalog.toml"


def test_catalog_builds_once_from_all_validated_tidy_publications(
    tmp_path: Path,
) -> None:
    resource_root = tmp_path / "resources"
    materialized = _snapshot(
        resource_root,
        database="chebi",
        asset="database",
        version="2026-07-07",
        publication=True,
    )
    direct_only = _snapshot(
        resource_root,
        database="eggnog",
        asset="mapper",
        version="5.0.2",
        publication=False,
    )
    manifest, output = _write_manifest(
        tmp_path,
        snapshots=[direct_only, materialized],
    )

    model = build_publication_catalog(manifest, output)
    first_bytes = output.read_bytes()
    repeated = build_publication_catalog(manifest, output)

    assert output.read_bytes() == first_bytes
    assert repeated == model
    parsed = tomllib.loads(first_bytes.decode("utf-8"))
    assert parsed["schema_version"] == "bioextract-publication-catalog-v1"
    assert parsed["resource_root"] == "../resources"
    assert parsed["summary"] == {
        "snapshot_count": 2,
        "publication_count": 1,
        "total_bytes": (resource_root / "chebi/database/2026-07-07/tidy/data.duckdb")
        .stat()
        .st_size,
    }
    assert len(parsed["publications"]) == 1
    publication = parsed["publications"][0]
    assert publication["database"] == "chebi"
    assert publication["snapshot_path"] == "chebi/database/2026-07-07"
    assert publication["path"] == "tidy/data.duckdb"
    assert publication["validation_status"] == "passed"
    assert publication["tables"] == [
        {"name": "mapping", "role": "canonical", "row_count": 1}
    ]


def test_catalog_cli_reports_build_summary(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    resource_root = tmp_path / "resources"
    snapshot = _snapshot(
        resource_root,
        database="rhea",
        asset="database",
        version="2026-06-10",
        publication=True,
    )
    manifest, output = _write_manifest(tmp_path, snapshots=[snapshot])

    main(["--manifest", str(manifest), "--output", str(output)])

    assert output.is_file()
    assert "snapshots=1 publications=1" in capsys.readouterr().out


def test_catalog_rejects_changed_snapshot_manifest(tmp_path: Path) -> None:
    resource_root = tmp_path / "resources"
    snapshot = _snapshot(
        resource_root,
        database="chebi",
        asset="database",
        version="2026-07-07",
        publication=True,
    )
    manifest, output = _write_manifest(tmp_path, snapshots=[snapshot])
    lock = resource_root / "chebi/database/2026-07-07/manifest.lock"
    lock.write_text("changed = true\n", encoding="utf-8")

    with pytest.raises(ValueError, match="snapshot manifest byte count"):
        build_publication_catalog(manifest, output)

    assert not output.exists()


def test_catalog_rejects_invalid_tidy_duckdb(tmp_path: Path) -> None:
    resource_root = tmp_path / "resources"
    snapshot = _snapshot(
        resource_root,
        database="chebi",
        asset="database",
        version="2026-07-07",
        publication=False,
    )
    tidy = resource_root / "chebi/database/2026-07-07/tidy"
    tidy.mkdir()
    with duckdb.connect(str(tidy / "data.duckdb")) as connection:
        connection.execute("CREATE TABLE data (id INTEGER)")
    manifest, output = _write_manifest(tmp_path, snapshots=[snapshot])

    with pytest.raises(ValueError, match="Invalid _bioextract table inventory"):
        build_publication_catalog(manifest, output)

    assert not output.exists()


def test_catalog_rejects_snapshot_path_traversal(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()
    manifest, output = _write_manifest(
        tmp_path,
        snapshots=[
            {
                "database": "bad",
                "asset": "database",
                "version": "1",
                "path": "../outside",
                "manifest": {
                    "path": "manifest.lock",
                    "sha256": "0" * 64,
                    "bytes": 0,
                },
            }
        ],
    )

    with pytest.raises(ValueError, match="canonical relative path"):
        build_publication_catalog(manifest, output)

    assert not output.exists()


def test_catalog_requires_biofetch_manifest_v1(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()
    meta = tmp_path / "meta"
    meta.mkdir()
    manifest = meta / "manifest.toml"
    manifest.write_text(
        "schema_version = 'unknown'\nresource_root = '../resources'\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="biofetch-manifest-v1"):
        build_publication_catalog(manifest, meta / "catalog.toml")
