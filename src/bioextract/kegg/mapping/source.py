from __future__ import annotations

import os
import re
import stat
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from bioextract.errors import CapabilityError, IntegrityError

from .constant import ORGANISM_ROLE_FILENAMES

_ORGANISM_CODE = re.compile(r"^[a-z]{3,4}$")
_RESERVED_DIRECTORIES = {"ko", "organism"}


@dataclass(frozen=True, slots=True)
class MappingSnapshot:
    """Retain one replayable KEGG mapping source or publication description."""

    mode: str
    source: Path | None = None
    organism_code: str | None = None
    organism_roles: Mapping[str, Path] | None = None
    organism_list: Path | None = None
    ko_pathway: Path | None = None
    organism_scope: tuple[str, ...] | None = None
    release_version: str | None = None
    publication_path: Path | None = None
    publication_capabilities: Mapping[str, bool] | None = None
    publication_members: tuple[str, ...] | None = None

    def with_organisms(self, organism_codes: Sequence[str]) -> MappingSnapshot:
        """Return an immutable source or publication scope."""
        codes = normalize_organism_scope(organism_codes)
        if self.mode == "files" and codes != (self.organism_code,):
            raise CapabilityError(
                "A single-organism KEGG mapping handle can only retain "
                f"{self.organism_code!r}"
            )
        if self.mode == "publication":
            available = set(self.publication_members or ())
            missing = sorted(set(codes) - available)
            if missing:
                raise CapabilityError(
                    f"KEGG mapping publication does not contain organisms: {missing}"
                )
        return replace(self, organism_scope=codes)


def create_directory_snapshot(
    source: os.PathLike[str] | str,
    *,
    organism_list: os.PathLike[str] | str | None,
    ko_pathway: os.PathLike[str] | str | None,
    release_version: str | None,
) -> MappingSnapshot:
    root = _require_directory(source, label="KEGG mapping source")
    organism_list_path = _resolve_optional_global(
        organism_list,
        discovered=root / "organism" / "list_organism.tsv",
        label="KEGG organism_list file",
    )
    ko_pathway_path = _resolve_optional_global(
        ko_pathway,
        discovered=root / "ko" / "ko_pathway.tsv",
        label="KEGG ko_pathway file",
    )
    _validate_distinct_paths(
        {
            "organism_list": organism_list_path,
            "ko_pathway": ko_pathway_path,
        }
    )
    return MappingSnapshot(
        mode="directory",
        source=root,
        organism_list=organism_list_path,
        ko_pathway=ko_pathway_path,
        release_version=_normalize_release_version(release_version),
    )


def create_files_snapshot(
    source: os.PathLike[str] | str | None,
    *,
    organism_code: str,
    gene_list: os.PathLike[str] | str | None,
    uniprot_conversion: os.PathLike[str] | str | None,
    ncbi_gene_conversion: os.PathLike[str] | str | None,
    gene_ko: os.PathLike[str] | str | None,
    gene_pathway: os.PathLike[str] | str | None,
    organism_list: os.PathLike[str] | str | None,
    ko_pathway: os.PathLike[str] | str | None,
    release_version: str | None,
) -> MappingSnapshot:
    code = validate_organism_code(organism_code)
    root = (
        None
        if source is None
        else _require_directory(source, label="KEGG organism mapping source")
    )
    supplied = {
        "gene_list": gene_list,
        "uniprot_conversion": uniprot_conversion,
        "ncbi_gene_conversion": ncbi_gene_conversion,
        "gene_ko": gene_ko,
        "gene_pathway": gene_pathway,
    }
    roles: dict[str, Path] = {}
    for role, filename in ORGANISM_ROLE_FILENAMES.items():
        value = supplied[role]
        if value is not None:
            roles[role] = _require_file(value, label=f"KEGG {role} file")
        elif root is not None and (candidate := root / filename).is_file():
            roles[role] = candidate
    if not roles:
        raise ValueError("KEGG mapping files require at least one organism role")
    organism_list_path = (
        None
        if organism_list is None
        else _require_file(organism_list, label="KEGG organism_list file")
    )
    ko_pathway_path = (
        None
        if ko_pathway is None
        else _require_file(ko_pathway, label="KEGG ko_pathway file")
    )
    _validate_distinct_paths(
        {**roles, "organism_list": organism_list_path, "ko_pathway": ko_pathway_path}
    )
    return MappingSnapshot(
        mode="files",
        source=root,
        organism_code=code,
        organism_roles=roles,
        organism_list=organism_list_path,
        ko_pathway=ko_pathway_path,
        release_version=_normalize_release_version(release_version),
    )


def normalize_organism_scope(organism_codes: Sequence[str]) -> tuple[str, ...]:
    codes = tuple(sorted({str(code).strip() for code in organism_codes}))
    if not codes:
        raise ValueError("organism_codes must contain at least one code")
    for code in codes:
        validate_organism_code(code)
    return codes


def validate_organism_code(organism_code: str) -> str:
    code = str(organism_code).strip()
    if _ORGANISM_CODE.fullmatch(code) is None:
        raise ValueError(
            "KEGG organism_code must match ^[a-z]{3,4}$ without case rewriting: "
            f"{organism_code!r}"
        )
    return code


def resolve_organism_work(
    snapshot: MappingSnapshot,
    *,
    validate_role_files: bool = True,
) -> tuple[tuple[str, Mapping[str, Path]], ...]:
    """Resolve one execution's organism work queue without recursive traversal."""
    if snapshot.mode == "files":
        code = snapshot.organism_code
        if code is None:
            raise IntegrityError("KEGG file snapshot is missing organism_code")
        return ((code, dict(snapshot.organism_roles or {})),)
    if snapshot.mode != "directory" or snapshot.source is None:
        raise CapabilityError("KEGG mapping source execution requires local files")
    if snapshot.organism_scope is not None:
        codes = snapshot.organism_scope
    else:
        codes_found: list[str] = []
        ignored: list[str] = []
        try:
            entries = tuple(os.scandir(snapshot.source))
        except OSError as error:
            raise IntegrityError(
                f"Cannot enumerate KEGG mapping source: {snapshot.source}"
            ) from error
        for entry in entries:
            if entry.name in _RESERVED_DIRECTORIES or not entry.is_dir(
                follow_symlinks=False
            ):
                continue
            if _ORGANISM_CODE.fullmatch(entry.name) is None:
                ignored.append(entry.name)
                continue
            codes_found.append(entry.name)
        if ignored:
            warnings.warn(
                "Ignored unsupported KEGG mapping directories: "
                + ", ".join(sorted(ignored)),
                stacklevel=2,
            )
        codes = tuple(sorted(set(codes_found)))
        if not codes:
            raise IntegrityError(
                f"KEGG mapping source contains no organism directories: {snapshot.source}"
            )
    work: list[tuple[str, Mapping[str, Path]]] = []
    for code in codes:
        directory = snapshot.source / code
        try:
            mode = directory.lstat().st_mode
        except OSError as error:
            raise IntegrityError(
                f"KEGG organism directory is unavailable: {directory}"
            ) from error
        if not stat.S_ISDIR(mode):
            raise IntegrityError(f"KEGG organism directory is unavailable: {directory}")
        roles: dict[str, Path] = {}
        for role, filename in ORGANISM_ROLE_FILENAMES.items():
            candidate = directory / filename
            if validate_role_files:
                try:
                    role_mode = candidate.lstat().st_mode
                except OSError as error:
                    raise IntegrityError(
                        f"KEGG organism {code!r} is missing required role "
                        f"{role!r}: {candidate}"
                    ) from error
                if not stat.S_ISREG(role_mode):
                    raise IntegrityError(
                        f"KEGG organism {code!r} role is not a regular file: "
                        f"{role!r}: {candidate}"
                    )
            roles[role] = candidate
        work.append((code, roles))
    return tuple(work)


def source_capabilities(snapshot: MappingSnapshot) -> dict[str, bool]:
    if snapshot.mode == "publication":
        return dict(snapshot.publication_capabilities or {})
    organism = (
        dict.fromkeys(ORGANISM_ROLE_FILENAMES, True)
        if snapshot.mode == "directory"
        else {
            role: role in (snapshot.organism_roles or {})
            for role in ORGANISM_ROLE_FILENAMES
        }
    )
    return {
        **organism,
        "organism_list": snapshot.organism_list is not None,
        "ko_pathway": snapshot.ko_pathway is not None,
    }


def _resolve_optional_global(
    explicit: os.PathLike[str] | str | None,
    *,
    discovered: Path,
    label: str,
) -> Path | None:
    if explicit is not None:
        return _require_file(explicit, label=label)
    return discovered if discovered.is_file() else None


def _require_file(value: os.PathLike[str] | str, *, label: str) -> Path:
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found or not a regular file: {path}")
    return path


def _require_directory(value: os.PathLike[str] | str, *, label: str) -> Path:
    path = Path(value)
    if not path.is_dir():
        raise FileNotFoundError(f"{label} not found or not a directory: {path}")
    return path


def _validate_distinct_paths(paths: Mapping[str, Path | None]) -> None:
    seen: dict[Path, str] = {}
    for role, path in paths.items():
        if path is None:
            continue
        resolved = path.resolve()
        if previous := seen.get(resolved):
            raise ValueError(
                f"KEGG source file is assigned to multiple roles: "
                f"{previous!r}, {role!r}: {path}"
            )
        seen[resolved] = role


def _normalize_release_version(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        raise ValueError("release_version must be non-empty when provided")
    return normalized
