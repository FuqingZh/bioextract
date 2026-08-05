from __future__ import annotations

import os
import tarfile
import tempfile
import zipfile
from collections.abc import Generator, Iterable, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import duckdb

from bioextract.errors import CapabilityError

from ._query import (
    RheaReactionSelection,
    _RheaPublication,  # pyright: ignore[reportPrivateUsage]  # sibling publication type
    create_group_selection,
    create_selection,
    open_rhea_publication,
)
from .constant import (
    COMPOUND_SOURCE_NAMES,
    CROSS_REFERENCE_SOURCE_NAMES,
    REACTION_SOURCE_NAMES,
    RELEASE_REQUIRED_SOURCES,
    SOURCE_BASENAMES,
    RheaNamespace,
)

__all__ = ["RheaDatabase"]


@dataclass(frozen=True, slots=True)
class RheaWriteResult:
    """Summary of one successfully committed Rhea DuckDB database.

    Examples:
        Inspect the immutable scope value:

        >>> report = RheaWriteResult(
        ...     path=Path("rhea.duckdb"),
        ...     scope="reactions",
        ...     tables=("reaction",),
        ...     row_counts={"reaction": 10},
        ...     source_files={"rdf": "rhea.rdf.gz"},
        ... )
        >>> report.scope
        'reactions'
    """

    path: Path
    scope: str
    tables: tuple[str, ...]
    row_counts: Mapping[str, int]
    source_files: Mapping[str, str]
    release_number: int | None = None
    release_date: str | None = None


@dataclass(frozen=True, slots=True)
class _RheaSnapshot:
    scope: str
    sources: Mapping[str, Path] = field(default_factory=dict[str, Path])
    release_source: Path | None = None


@dataclass(slots=True)
class RheaDatabase:
    """Build a query-ready DuckDB database from local Rhea release files.

    :meth:`from_files` exposes independently useful or mixed reaction,
    compound, and cross-reference capabilities. :meth:`from_release` accepts
    either an extracted release directory or a ``zip``/``tar`` archive and
    requires the complete official release asset set.

    Examples:
        Open only the reaction component:

        >>> db = RheaDatabase.from_files(
        ...     rdf="rhea.rdf.gz",
        ...     directions="rhea-directions.tsv",
        ... )
        >>> db.snapshot.scope
        'reactions'
    """

    snapshot: _RheaSnapshot
    _publication: _RheaPublication | None = field(default=None, repr=False)

    @classmethod
    def from_files(
        cls,
        *,
        rdf: os.PathLike[str] | str | None = None,
        directions: os.PathLike[str] | str | None = None,
        relationships: os.PathLike[str] | str | None = None,
        obsolete_reactions: os.PathLike[str] | str | None = None,
        reaction_smiles: os.PathLike[str] | str | None = None,
        sdf: os.PathLike[str] | str | None = None,
        chebi_names: os.PathLike[str] | str | None = None,
        chebi_ph7_3_mapping: os.PathLike[str] | str | None = None,
        xrefs: os.PathLike[str] | str | None = None,
        uniprot_sprot: os.PathLike[str] | str | None = None,
        uniprot_trembl: os.PathLike[str] | str | None = None,
    ) -> RheaDatabase:
        """Create a capability-scoped handle from explicit Rhea files.

        A reaction role requires both ``rdf`` and ``directions``. Compound and
        cross-reference roles are independently constructible. A handle with
        exactly one role group has that group's scope; mixed groups have
        ``partial`` scope.

        Args:
            rdf: Optional official ``rhea.rdf`` or ``rhea.rdf.gz``.
            directions: Optional official reaction direction quartet table.
            relationships: Optional reaction hierarchy table.
            obsolete_reactions: Optional obsolete reaction ID table.
            reaction_smiles: Optional headerless reaction SMILES table.
            sdf: Optional Rhea compound structure SDF.
            chebi_names: Optional headerless ChEBI ID/name table.
            chebi_ph7_3_mapping: Optional pH 7.3 ChEBI mapping table.
            xrefs: Optional aggregate Rhea cross-reference table.
            uniprot_sprot: Optional reviewed UniProt reaction mapping.
            uniprot_trembl: Optional unreviewed UniProt mapping, plain or
                gzip-compressed.

        Returns:
            A source-backed handle containing only the supplied capabilities.

        Raises:
            FileNotFoundError: If a supplied path does not exist.
            ValueError: If no role is supplied, a reaction dependency is
                missing, a path is not a file, or two roles resolve to the same
                physical file.

        Notes:
            ``obsolete_reactions`` is the caller-facing role. It intentionally
            maps to the stable internal and provenance role ``obsoletes``.

        Examples:
            Compound inputs remain independently constructible:

            >>> db = RheaDatabase.from_files(chebi_names="names.tsv")  # doctest: +SKIP
            >>> db.snapshot.scope  # doctest: +SKIP
            'compounds'
        """
        explicit_values = {
            "rdf": rdf,
            "directions": directions,
            "relationships": relationships,
            "obsolete_reactions": obsolete_reactions,
            "reaction_smiles": reaction_smiles,
            "sdf": sdf,
            "chebi_names": chebi_names,
            "chebi_ph7_3_mapping": chebi_ph7_3_mapping,
            "xrefs": xrefs,
            "uniprot_sprot": uniprot_sprot,
            "uniprot_trembl": uniprot_trembl,
        }
        supplied = {
            name for name, value in explicit_values.items() if value is not None
        }
        if not supplied:
            raise ValueError("At least one Rhea input file role must be provided")

        reaction_roles = {
            "rdf",
            "directions",
            "relationships",
            "obsolete_reactions",
            "reaction_smiles",
        }
        if supplied & reaction_roles and not {"rdf", "directions"} <= supplied:
            missing = sorted({"rdf", "directions"} - supplied)
            raise ValueError(
                "Rhea reaction input roles require both rdf and directions; "
                f"missing roles: {', '.join(missing)}"
            )

        groups = {
            "reactions": reaction_roles,
            "compounds": {"sdf", "chebi_names", "chebi_ph7_3_mapping"},
            "cross_references": {"xrefs", "uniprot_sprot", "uniprot_trembl"},
        }
        present_groups = [name for name, roles in groups.items() if supplied & roles]
        scope = present_groups[0] if len(present_groups) == 1 else "partial"
        sources = _validate_explicit_sources(explicit_values)
        if "obsolete_reactions" in sources:
            sources["obsoletes"] = sources.pop("obsolete_reactions")
        return cls(snapshot=_RheaSnapshot(scope=scope, sources=sources))

    @classmethod
    def from_release(
        cls,
        source: os.PathLike[str] | str,
    ) -> RheaDatabase:
        """Create a strict handle from a complete Rhea release.

        ``source`` may be an extracted directory, a zip archive, or a tar
        archive. Compression is detected from content rather than filename.

        Examples:
            A missing release path is rejected immediately:

            >>> try:
            ...     RheaDatabase.from_release("missing-rhea-release")
            ... except FileNotFoundError as error:
            ...     print(error.filename)
            missing-rhea-release
        """
        path = _validate_path(source)
        if path.is_dir():
            sources = _discover_release_files(path)
            _validate_complete_release(sources)
            return cls(snapshot=_RheaSnapshot(scope="release", sources=sources))

        _inspect_archive_release(path)
        return cls(snapshot=_RheaSnapshot(scope="release", release_source=path))

    @classmethod
    def from_duckdb(
        cls,
        path: os.PathLike[str] | str,
    ) -> RheaDatabase:
        """Open a validated bioextract Rhea publication for domain queries.

        The returned handle opens short-lived read-only DuckDB connections for
        extraction. It cannot be passed back to :meth:`write_duckdb`.

        Args:
            path: Existing bioextract Rhea DuckDB publication.

        Returns:
            A publication-backed handle for reaction selections.

        Raises:
            FileNotFoundError: If ``path`` does not exist.
            ValueError: If the file is not a supported bioextract Rhea
                publication or its recorded inventory is inconsistent.

        Examples:
            Open a previously published database for read-only selection:

            >>> db = RheaDatabase.from_duckdb("rhea.duckdb")  # doctest: +SKIP
            >>> db.snapshot.scope  # doctest: +SKIP
            'publication'
        """
        publication = open_rhea_publication(Path(path))
        return cls(
            snapshot=_RheaSnapshot(scope="publication"),
            _publication=publication,
        )

    def select_reactions(
        self,
        ids: Iterable[str],
        *,
        namespace: RheaNamespace,
        include_obsolete: bool = False,
    ) -> RheaReactionSelection:
        """Create a deferred exact-reaction selection from one ID namespace.

        Args:
            ids: Identifiers from the declared namespace.
            namespace: Rhea, ChEBI, UniProt, EC, GO, or another supported
                official Rhea cross-reference namespace.
            include_obsolete: Include obsolete reaction records when available.

        Returns:
            An immutable query plan with eager ``extract_*`` terminals.

        Raises:
            CapabilityError: If this is not a publication-backed handle or
                its partial publication lacks the required relations.
            ValueError: If ``namespace`` or an identifier is invalid.

        Examples:
            Defer a ChEBI-to-reaction lookup until extraction:

            >>> selection = db.select_reactions(  # doctest: +SKIP
            ...     ["CHEBI:15377"],
            ...     namespace="chebi",
            ... )
            >>> selection.namespace  # doctest: +SKIP
            'chebi'
        """
        return create_selection(
            self,
            ids,
            namespace=namespace,
            include_obsolete=include_obsolete,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: RheaNamespace,
        include_obsolete: bool = False,
    ) -> RheaReactionSelection:
        """Create a deferred reaction selection preserving group isolation.

        Args:
            ids_by_group: Mapping of non-empty group IDs to identifiers.
            namespace: Namespace shared by every identifier in this selection.
            include_obsolete: Include obsolete reaction records when available.

        Returns:
            An immutable grouped query plan with ``group_id`` lineage.

        Examples:
            Retain group identity while selecting by EC number:

            >>> selection = db.select_groups(  # doctest: +SKIP
            ...     {"glycolysis": ["1.2.1.12"]},
            ...     namespace="ec",
            ... )
            >>> selection.extract_matches().columns[:2]  # doctest: +SKIP
            ['group_id', 'input_id']
        """
        return create_group_selection(
            self,
            ids_by_group,
            namespace=namespace,
            include_obsolete=include_obsolete,
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: Literal["fail", "replace"] = "fail",
        include_source_hashes: bool = False,
    ) -> RheaWriteResult:
        """Atomically write the selected Rhea components to one DuckDB file.

        Args:
            path: Destination database path.
            if_exists: ``"fail"`` preserves an existing destination;
                ``"replace"`` atomically replaces it after a successful build.
            include_source_hashes: Calculate SHA-256 values for source provenance.

        Returns:
            The committed database path, table inventory, counts, and release
            metadata.

        Examples:
            Invalid overwrite policy values are rejected before writing:

            >>> db = object.__new__(RheaDatabase)
            >>> try:
            ...     db.write_duckdb("rhea.duckdb", if_exists="skip")
            ... except ValueError as error:
            ...     print(str(error))
            if_exists must be 'fail' or 'replace'
        """
        if if_exists not in {"fail", "replace"}:
            raise ValueError("if_exists must be 'fail' or 'replace'")
        if self._publication is not None:
            raise RuntimeError(
                "A publication-backed RheaDatabase cannot be republished; "
                "construct a handle from official source files"
            )

        from ._duckdb import write_rhea_duckdb

        destination = Path(path)
        with self._resolved_sources() as (sources, display_paths):
            return write_rhea_duckdb(
                sources=sources,
                display_paths=display_paths,
                scope=self.snapshot.scope,
                path=destination,
                if_exists=if_exists,
                include_source_hashes=include_source_hashes,
            )

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Open a new native read-only connection to this publication.

        The caller owns the returned connection and should close it or use it
        as a context manager. Each call returns an independent connection.

        Examples:
            Inspect a publication through native read-only SQL:

            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.execute(
            ...         "SELECT count(*) FROM reaction"
            ...     ).fetchone()[0]
            >>> count > 0  # doctest: +SKIP
            True
        """
        publication = self._require_publication()
        return duckdb.connect(str(publication.path), read_only=True)

    def _require_publication(self) -> _RheaPublication:
        if self._publication is None:
            raise CapabilityError(
                "Rhea domain selection requires a publication-backed handle; "
                "write a DuckDB publication and open it with "
                "RheaDatabase.from_duckdb()"
            )
        return self._publication

    @contextmanager
    def _resolved_sources(
        self,
    ) -> Generator[tuple[Mapping[str, Path], Mapping[str, str]]]:
        if self.snapshot.release_source is None:
            yield (
                self.snapshot.sources,
                {name: str(path) for name, path in self.snapshot.sources.items()},
            )
            return

        archive = self.snapshot.release_source
        with tempfile.TemporaryDirectory(prefix="bioextract-rhea-") as dir_tmp:
            root = Path(dir_tmp)
            _extract_archive(archive, root)
            sources = _discover_release_files(root)
            _validate_complete_release(sources)
            display_paths = {
                name: f"{archive}::{path.relative_to(root)}"
                for name, path in sources.items()
            }
            yield sources, display_paths


def _validate_explicit_sources(
    values: Mapping[str, os.PathLike[str] | str | None],
) -> dict[str, Path]:
    sources: dict[str, Path] = {}
    physical_roles: dict[tuple[int, int], tuple[str, Path]] = {}
    for name, value in values.items():
        if value is None:
            continue
        path = _validate_path(value).resolve()
        if not path.is_file():
            raise ValueError(f"Rhea input is not a file: {path}")
        stat = path.stat()
        physical_id = (stat.st_dev, stat.st_ino)
        if physical_id in physical_roles:
            other_name, other_path = physical_roles[physical_id]
            raise ValueError(
                f"Rhea input roles {other_name!r} and {name!r} refer to the "
                f"same physical file: {other_path}"
            )
        physical_roles[physical_id] = (name, path)
        sources[name] = path
    return sources


def _validate_path(value: os.PathLike[str] | str) -> Path:
    path = Path(value)
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def _discover_release_files(root: Path) -> dict[str, Path]:
    candidates_by_basename: dict[str, list[Path]] = {}
    for path in root.rglob("*"):
        if path.is_file():
            candidates_by_basename.setdefault(path.name, []).append(path)

    sources: dict[str, Path] = {}
    for name, basenames in SOURCE_BASENAMES.items():
        matches = [
            path
            for basename in basenames
            for path in candidates_by_basename.get(basename, [])
        ]
        if len(matches) > 1:
            raise ValueError(
                f"Rhea release contains multiple candidates for {name}: "
                + ", ".join(str(path) for path in matches)
            )
        if matches:
            sources[name] = matches[0]
    return sources


def _validate_complete_release(
    sources: Mapping[str, Path],
) -> None:
    missing = [name for name in RELEASE_REQUIRED_SOURCES if name not in sources]
    if missing:
        expected = {name: SOURCE_BASENAMES[name][0] for name in missing}
        raise ValueError(f"Incomplete Rhea release; missing sources: {expected}")


def _inspect_archive_release(path: Path) -> None:
    names = _archive_member_names(path)
    basenames = [Path(name).name for name in names if not name.endswith("/")]
    missing: list[str] = []
    for name in RELEASE_REQUIRED_SOURCES:
        matches = [
            basename for basename in basenames if basename in SOURCE_BASENAMES[name]
        ]
        if len(matches) > 1:
            raise ValueError(
                f"Rhea release archive contains multiple candidates for "
                f"{name}: {matches}"
            )
        if not matches:
            missing.append(name)
    if missing:
        raise ValueError(f"Incomplete Rhea release archive; missing sources: {missing}")


def _archive_member_names(path: Path) -> list[str]:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            return archive.namelist()
    if tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            return [member.name for member in archive.getmembers()]
    raise ValueError(f"Unsupported Rhea release archive: {path}")


def _extract_archive(path: Path, root: Path) -> None:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for member in archive.infolist():
                target = (root / member.filename).resolve()
                if not target.is_relative_to(root.resolve()):
                    raise ValueError(
                        f"Unsafe path in Rhea release archive: {member.filename}"
                    )
            archive.extractall(root)
        return
    if tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            archive.extractall(root, filter="data")
        return
    raise ValueError(f"Unsupported Rhea release archive: {path}")


assert set(REACTION_SOURCE_NAMES).isdisjoint(COMPOUND_SOURCE_NAMES)
assert set(COMPOUND_SOURCE_NAMES).isdisjoint(CROSS_REFERENCE_SOURCE_NAMES)
