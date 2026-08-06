from __future__ import annotations

import gzip
import io
import os
import shutil
import tarfile
import tempfile
import zipfile
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Literal, TextIO, cast

import duckdb
import polars as pl

from bioextract._publication import (
    DuckDBWriteResult,
    RelationSpec,
    SourceFileRecord,
    write_duckdb_publication,
)
from bioextract.errors import CapabilityError

from ._canonical import build_canonical_relations
from ._query import (
    SOURCE_SCHEMA_PROFILE,
    ChEBICompoundSelection,
    _ChEBIPublication,  # pyright: ignore[reportPrivateUsage]  # sibling publication type
    create_group_selection,
    create_selection,
    open_chebi_publication,
)

__all__ = ["ChEBIDatabase"]

_SCHEMA_VERSION = "chebi-duckdb-v1"
_TABLE_FILES = {
    "compound": "compounds.tsv",
    "compound_name": "names.tsv",
    "compound_relation": "relation.tsv",
    "secondary_id": "secondary_ids.tsv",
    "database_accession": "database_accession.tsv",
    "structure": "structures.tsv",
    "chemical_data": "chemical_data.tsv",
}
_TABLE_ROLES = {
    "compound": "entity",
    "compound_name": "annotation",
    "compound_relation": "relationship",
    "secondary_id": "identifier_mapping",
    "database_accession": "cross_reference",
    "structure": "structure",
    "chemical_data": "property",
}


@dataclass(frozen=True, slots=True)
class _ChEBISnapshot:
    table_sources: Mapping[str, Path] = field(default_factory=dict[str, Path])
    release_source: Path | None = None
    file_obo: Path | None = None
    file_sdf: Path | None = None
    file_chemont_obo: Path | None = None


@dataclass(slots=True)
class ChEBIDatabase:
    """Build one relational DuckDB from local ChEBI and ChemOnt snapshots.

    Official ChEBI TSV columns retain their source headers. OBO-derived ChEBI
    and ChemOnt relations use stable snake_case columns and remain separate
    relations inside the same database.

    Examples:
        Construct an ontology-only handle from a local ChEBI OBO snapshot:

        >>> db = ChEBIDatabase.from_obo("data/chebi.obo")
        >>> db.snapshot.file_obo.name
        'chebi.obo'
    """

    snapshot: _ChEBISnapshot
    _publication: _ChEBIPublication | None = field(default=None, repr=False)

    @classmethod
    def from_table_files(
        cls,
        source: os.PathLike[str] | str | None = None,
        *,
        compounds: os.PathLike[str] | str | None = None,
        names: os.PathLike[str] | str | None = None,
        relations: os.PathLike[str] | str | None = None,
        secondary_ids: os.PathLike[str] | str | None = None,
        database_accessions: os.PathLike[str] | str | None = None,
        structures: os.PathLike[str] | str | None = None,
        chemical_data: os.PathLike[str] | str | None = None,
        chemont_obo: os.PathLike[str] | str | None = None,
    ) -> ChEBIDatabase:
        """Open explicit ChEBI table files for a partial or complete build.

        Each parameter names the logical role of an official ChEBI file.
        Plain and gzip-compressed TSV inputs are detected by content.

        Examples:
            Missing required compounds input fails at construction:

            >>> try:
            ...     ChEBIDatabase.from_table_files(
            ...         compounds="missing-compounds.tsv"
            ...     )
            ... except FileNotFoundError as error:
            ...     print(error.filename)
            missing-compounds.tsv
        """
        values = {
            "compound": compounds,
            "compound_name": names,
            "compound_relation": relations,
            "secondary_id": secondary_ids,
            "database_accession": database_accessions,
            "structure": structures,
            "chemical_data": chemical_data,
        }
        explicit_sources = {
            table_name: _require_file(path)
            for table_name, path in values.items()
            if path is not None
        }
        _validate_table_source_inventory(explicit_sources)
        chemont = _optional_file(chemont_obo)
        if source is None:
            _require_compounds(explicit_sources)
            return cls(
                snapshot=_ChEBISnapshot(
                    table_sources=explicit_sources,
                    file_chemont_obo=chemont,
                )
            )

        source_path = _require_file_or_directory(source)
        if source_path.is_dir():
            discovered = _discover_table_files(
                source_path,
                skipped_roles=explicit_sources,
            )
            table_sources = {**discovered, **explicit_sources}
            _require_compounds(table_sources)
            _validate_table_source_inventory(table_sources)
            return cls(
                snapshot=_ChEBISnapshot(
                    table_sources=table_sources,
                    file_chemont_obo=chemont,
                )
            )

        _inspect_release_archive(
            source_path,
            skipped_roles=explicit_sources,
        )
        return cls(
            snapshot=_ChEBISnapshot(
                table_sources=explicit_sources,
                release_source=source_path,
                file_chemont_obo=chemont,
            )
        )

    @classmethod
    def from_obo(
        cls,
        source: os.PathLike[str] | str,
        *,
        sdf: os.PathLike[str] | str | None = None,
        chemont_obo: os.PathLike[str] | str | None = None,
    ) -> ChEBIDatabase:
        """Open ChEBI OBO and optional ChemOnt OBO input.

        Ordinary files, gzip streams, zip archives, and tar archives are
        recognized from their content; compression suffixes are not required.

        Examples:
            A missing ontology input fails at construction:

            >>> try:
            ...     ChEBIDatabase.from_obo("missing-chebi.obo")
            ... except FileNotFoundError as error:
            ...     print(error.filename)
            missing-chebi.obo
        """
        source_path = _require_file_or_directory(source)
        if source_path.is_dir():
            file_obo = _discover_ontology_file(source_path)
            _inspect_ontology_source(file_obo)
            file_sdf = (
                _optional_file(sdf)
                if sdf is not None
                else _discover_supplement_file(source_path, suffix=".sdf")
            )
        else:
            _inspect_ontology_source(source_path)
            file_obo = source_path
            if sdf is not None:
                file_sdf = _optional_file(sdf)
            elif zipfile.is_zipfile(source_path) or tarfile.is_tarfile(source_path):
                file_sdf = (
                    source_path
                    if _archive_has_exactly_one(source_path, suffix=".sdf")
                    else None
                )
            else:
                file_sdf = None
        chemont = _optional_file(chemont_obo)
        if chemont is not None:
            _inspect_ontology_source(chemont)
        return cls(
            snapshot=_ChEBISnapshot(
                file_obo=file_obo,
                file_sdf=file_sdf,
                file_chemont_obo=chemont,
            )
        )

    @classmethod
    def from_duckdb(
        cls,
        path: os.PathLike[str] | str,
    ) -> ChEBIDatabase:
        """Open a validated ChEBI publication for read-only domain access.

        Examples:
            Open a publication before selecting canonical compounds:

            >>> db = ChEBIDatabase.from_duckdb(  # doctest: +SKIP
            ...     "chebi.duckdb"
            ... )
            >>> type(db).__name__  # doctest: +SKIP
            'ChEBIDatabase'
        """
        publication = open_chebi_publication(Path(path))
        return cls(
            snapshot=_ChEBISnapshot(),
            _publication=publication,
        )

    def connect(self) -> duckdb.DuckDBPyConnection:
        """Return a new native read-only DuckDB connection.

        The caller owns the connection and may execute arbitrary read-only SQL
        against both domain relations and ``_bioextract`` provenance.

        Examples:
            Run native SQL while retaining publication immutability:

            >>> with db.connect() as connection:  # doctest: +SKIP
            ...     count = connection.execute(
            ...         "SELECT count(*) FROM compound"
            ...     ).fetchone()[0]
            >>> count > 0  # doctest: +SKIP
            True
        """
        publication = self._require_publication()
        return duckdb.connect(str(publication.path), read_only=True)

    def select_compounds(
        self,
        ids: Iterable[str],
        *,
        namespace: str,
        min_star_rating: int = 1,
        include_obsolete: bool = False,
    ) -> ChEBICompoundSelection:
        """Create a deferred exact compound selection.

        Examples:
            Resolve primary and secondary ChEBI identifiers together:

            >>> selection = db.select_compounds(  # doctest: +SKIP
            ...     ["CHEBI:15377", "CHEBI:10743"],
            ...     namespace="chebi",
            ... )
            >>> selection.extract_matches()["chebi_id"].unique().to_list()  # doctest: +SKIP
            ['CHEBI:15377']
        """
        return create_selection(
            self,
            ids,
            namespace=namespace,
            min_star_rating=min_star_rating,
            include_obsolete=include_obsolete,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
        *,
        namespace: str,
        min_star_rating: int = 1,
        include_obsolete: bool = False,
    ) -> ChEBICompoundSelection:
        """Create a deferred grouped compound selection.

        Examples:
            Preserve group lineage across compound resolution:

            >>> selection = db.select_groups(  # doctest: +SKIP
            ...     {"solvent": ["CHEBI:15377"]},
            ...     namespace="chebi",
            ... )
            >>> selection.extract_matches()["group_id"].to_list()  # doctest: +SKIP
            ['solvent']
        """
        return create_group_selection(
            self,
            ids_by_group,
            namespace=namespace,
            min_star_rating=min_star_rating,
            include_obsolete=include_obsolete,
        )

    def write_duckdb(
        self,
        path: os.PathLike[str] | str,
        *,
        if_exists: Literal["fail", "replace"] = "fail",
    ) -> DuckDBWriteResult:
        """Atomically publish the selected ChEBI and ChemOnt relations.

        Args:
            path: Destination DuckDB file.
            if_exists: Preserve or atomically replace an existing destination.

        Examples:
            Invalid overwrite policies fail before publication:

            >>> db = object.__new__(ChEBIDatabase)
            >>> try:
            ...     db.write_duckdb("chebi.duckdb", if_exists="skip")
            ... except ValueError as error:
            ...     print(str(error))
            if_exists must be 'fail' or 'replace'
        """
        if if_exists not in {"fail", "replace"}:
            raise ValueError("if_exists must be 'fail' or 'replace'")
        if self._publication is not None:
            raise RuntimeError(
                "A publication-backed ChEBIDatabase cannot be republished; "
                "construct a handle from official source files"
            )

        with _prepared_table_sources(self.snapshot) as (
            table_sources,
            provenance_sources,
        ):
            relations: list[RelationSpec] = []
            validation_issues = ()
            if self.snapshot.file_obo is not None:
                canonical = build_canonical_relations(
                    self.snapshot.file_obo,
                    file_sdf=self.snapshot.file_sdf,
                )
                relations.extend(canonical.relations)
                validation_issues = canonical.validation_issues
                provenance_sources.append(
                    _source_record("chebi_obo", self.snapshot.file_obo)
                )
                if self.snapshot.file_sdf is not None:
                    provenance_sources.append(
                        _source_record("chebi_sdf", self.snapshot.file_sdf)
                    )
            else:
                relations.extend(
                    RelationSpec(
                        table_name=table_name,
                        frame=(frame := _scan_official_tsv(file_source)),
                        role=_TABLE_ROLES[table_name],
                        source_columns=tuple(frame.collect_schema().names()),
                    )
                    for table_name, file_source in table_sources.items()
                )
            if self.snapshot.file_chemont_obo is not None:
                relations.extend(
                    _obo_relations(
                        self.snapshot.file_chemont_obo,
                        namespace="chemont",
                    )
                )
                provenance_sources.append(
                    _source_record(
                        "chemont_obo",
                        self.snapshot.file_chemont_obo,
                    )
                )
            if not relations:
                raise ValueError("No ChEBI or ChemOnt relations were selected")

            return write_duckdb_publication(
                relations,
                path,
                resource_name="chebi",
                resource_schema_version=_SCHEMA_VERSION,
                source_schema_profile=SOURCE_SCHEMA_PROFILE,
                sources=provenance_sources,
                scope=_scope(self.snapshot, bool(table_sources)),
                release_version=None,
                if_exists=if_exists,
                validation_issues=validation_issues,
            )

    def _require_publication(self) -> _ChEBIPublication:
        if self._publication is None:
            raise CapabilityError(
                "ChEBI domain selection requires a publication-backed handle; "
                "write a DuckDB publication and open it with "
                "ChEBIDatabase.from_duckdb()"
            )
        return self._publication


def _scan_official_tsv(file_source: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file_source,
        separator="\t",
        has_header=True,
        infer_schema_length=10_000,
        null_values=[""],
        quote_char='"',
    )


def _obo_relations(file_obo: Path, *, namespace: str) -> list[RelationSpec]:
    frames = _read_obo_frames(file_obo)
    return [
        RelationSpec(
            table_name=f"{namespace}_{relation_name}",
            frame=frame.lazy(),
            role=role,
        )
        for relation_name, frame, role in (
            ("term", frames["term"], "entity"),
            ("term_relation", frames["term_relation"], "relationship"),
            ("term_synonym", frames["term_synonym"], "annotation"),
            ("term_xref", frames["term_xref"], "cross_reference"),
        )
    ]


def _read_obo_frames(file_obo: Path) -> dict[str, pl.DataFrame]:
    terms: list[dict[str, object]] = []
    relations: list[dict[str, str]] = []
    synonyms: list[dict[str, str]] = []
    xrefs: list[dict[str, str]] = []
    stanza: dict[str, list[str]] | None = None

    def flush() -> None:
        if stanza is None or "id" not in stanza:
            return
        term_id = stanza["id"][0]
        terms.append(
            {
                "term_id": term_id,
                "term_name": _first(stanza, "name"),
                "definition": _quoted(_first(stanza, "def")),
                "is_obsolete": _first(stanza, "is_obsolete") == "true",
            }
        )
        for parent in stanza.get("is_a", []):
            relations.append(
                {
                    "child_term_id": term_id,
                    "parent_term_id": parent.split()[0],
                    "relation_type": "is_a",
                }
            )
        for relationship in stanza.get("relationship", []):
            fields = relationship.split()
            if len(fields) >= 2:
                relations.append(
                    {
                        "child_term_id": term_id,
                        "parent_term_id": fields[1],
                        "relation_type": fields[0],
                    }
                )
        for synonym in stanza.get("synonym", []):
            synonyms.append(
                {
                    "term_id": term_id,
                    "synonym": _quoted(synonym) or "",
                    "scope": _synonym_scope(synonym),
                }
            )
        for xref in stanza.get("xref", []):
            xrefs.append(
                {
                    "term_id": term_id,
                    "xref_id": xref.split()[0],
                }
            )

    with _open_ontology_text(file_obo) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "[Term]":
                flush()
                stanza = {}
                continue
            if line.startswith("["):
                flush()
                stanza = None
                continue
            if stanza is None or not line or line.startswith("!"):
                continue
            key, separator, value = line.partition(":")
            if separator:
                stanza.setdefault(key, []).append(value.strip())
        flush()

    return {
        "term": pl.DataFrame(
            terms,
            schema={
                "term_id": pl.String,
                "term_name": pl.String,
                "definition": pl.String,
                "is_obsolete": pl.Boolean,
            },
        ),
        "term_relation": pl.DataFrame(
            relations,
            schema={
                "child_term_id": pl.String,
                "parent_term_id": pl.String,
                "relation_type": pl.String,
            },
        ),
        "term_synonym": pl.DataFrame(
            synonyms,
            schema={
                "term_id": pl.String,
                "synonym": pl.String,
                "scope": pl.String,
            },
        ),
        "term_xref": pl.DataFrame(
            xrefs,
            schema={
                "term_id": pl.String,
                "xref_id": pl.String,
            },
        ),
    }


@contextmanager
def _prepared_table_sources(
    snapshot: _ChEBISnapshot,
) -> Generator[tuple[dict[str, Path], list[SourceFileRecord]]]:
    with tempfile.TemporaryDirectory(prefix="bioextract-chebi-") as dir_tmp_value:
        dir_tmp = Path(dir_tmp_value)
        if snapshot.release_source is not None:
            discovered = _extract_release_archive(
                snapshot.release_source,
                dir_tmp,
                skipped_roles=snapshot.table_sources,
            )
            table_sources = {**discovered, **snapshot.table_sources}
            _require_compounds(table_sources)
            _validate_table_source_inventory(table_sources)
            provenance = [_source_record("release_archive", snapshot.release_source)]
            provenance.extend(
                _source_record(table_name, file_source)
                for table_name, file_source in snapshot.table_sources.items()
            )
        else:
            table_sources = dict(snapshot.table_sources)
            provenance = [
                _source_record(table_name, file_source)
                for table_name, file_source in table_sources.items()
            ]

        prepared: dict[str, Path] = {}
        for table_name, file_source in table_sources.items():
            if _is_gzip(file_source):
                file_plain = dir_tmp / f"{table_name}.tsv"
                with (
                    gzip.open(file_source, "rb") as handle_in,
                    file_plain.open("wb") as handle_out,
                ):
                    shutil.copyfileobj(handle_in, handle_out)
                prepared[table_name] = file_plain
            else:
                prepared[table_name] = file_source
        yield prepared, provenance


def _discover_table_files(
    directory: Path,
    *,
    skipped_roles: Mapping[str, Path] | Iterable[str] = (),
) -> dict[str, Path]:
    skipped = set(skipped_roles)
    by_base_name: dict[str, list[Path]] = {}
    for candidate in directory.rglob("*"):
        if candidate.is_file():
            name = candidate.name.removesuffix(".gz").lower()
            by_base_name.setdefault(name, []).append(candidate)
    discovered: dict[str, Path] = {}
    for table_name, base_name in _TABLE_FILES.items():
        if table_name in skipped:
            continue
        matches = by_base_name.get(base_name.lower(), [])
        if len(matches) > 1:
            raise ValueError(
                f"ChEBI source contains multiple candidates for {table_name}: "
                + ", ".join(str(path) for path in sorted(matches))
            )
        if matches:
            discovered[table_name] = matches[0]
    return discovered


def _discover_ontology_file(directory: Path) -> Path:
    matches = sorted(
        candidate
        for candidate in directory.rglob("*")
        if candidate.is_file()
        and candidate.name.lower().removesuffix(".gz").endswith(".obo")
    )
    if len(matches) != 1:
        raise ValueError(
            "ChEBI source directory must contain exactly one .obo candidate; "
            f"found {len(matches)}"
        )
    return matches[0]


def _discover_supplement_file(directory: Path, *, suffix: str) -> Path | None:
    matches = sorted(
        candidate
        for candidate in directory.rglob("*")
        if candidate.is_file()
        and candidate.name.lower().removesuffix(".gz").endswith(suffix)
    )
    if len(matches) > 1:
        raise ValueError(
            f"ChEBI source directory contains multiple {suffix} candidates: "
            + ", ".join(str(path) for path in matches)
        )
    return matches[0] if matches else None


def _archive_has_exactly_one(path: Path, *, suffix: str) -> bool:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            names = [
                member.filename for member in archive.infolist() if not member.is_dir()
            ]
    elif tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            names = [member.name for member in archive.getmembers() if member.isfile()]
    else:
        return False
    for member in names:
        _validate_archive_member(member)
    matches = [
        member
        for member in names
        if member.lower().removesuffix(".gz").endswith(suffix)
    ]
    if len(matches) > 1:
        raise ValueError(
            f"ChEBI source archive must contain at most one {suffix} member; "
            f"found {len(matches)}"
        )
    return bool(matches)


def _inspect_release_archive(
    path: Path,
    *,
    skipped_roles: Mapping[str, Path] | Iterable[str] = (),
) -> None:
    with tempfile.TemporaryDirectory(prefix="bioextract-chebi-inspect-") as value:
        extracted = _extract_release_archive(
            path,
            Path(value),
            skipped_roles=skipped_roles,
        )
        if "compound" not in extracted and "compound" not in set(skipped_roles):
            _require_compounds(extracted)


def _extract_release_archive(
    path: Path,
    directory: Path,
    *,
    skipped_roles: Mapping[str, Path] | Iterable[str] = (),
) -> dict[str, Path]:
    expected_by_basename = {
        base_name.lower(): table_name for table_name, base_name in _TABLE_FILES.items()
    }
    expected_by_basename.update(
        {
            f"{base_name.lower()}.gz": table_name
            for table_name, base_name in _TABLE_FILES.items()
        }
    )
    extracted: dict[str, Path] = {}
    skipped = set(skipped_roles)
    if zipfile.is_zipfile(path):
        _extract_zip_release(
            path,
            expected_by_basename,
            directory,
            extracted,
            skipped,
        )
    elif tarfile.is_tarfile(path):
        _extract_tar_release(
            path,
            expected_by_basename,
            directory,
            extracted,
            skipped,
        )
    else:
        raise ValueError(
            f"ChEBI release source is not a directory, zip, or tar archive: {path}"
        )
    return extracted


def _extract_zip_release(
    path: Path,
    expected_by_basename: Mapping[str, str],
    directory: Path,
    extracted: dict[str, Path],
    skipped_roles: set[str],
) -> None:
    with zipfile.ZipFile(path) as archive:
        members = (
            (
                member.filename,
                cast(
                    Callable[[], BinaryIO],
                    lambda member=member: archive.open(member.filename),
                ),
            )
            for member in archive.infolist()
            if not member.is_dir()
        )
        _copy_release_members(
            members,
            expected_by_basename,
            directory,
            extracted,
            skipped_roles,
        )


def _extract_tar_release(
    path: Path,
    expected_by_basename: Mapping[str, str],
    directory: Path,
    extracted: dict[str, Path],
    skipped_roles: set[str],
) -> None:
    with tarfile.open(path, mode="r:*") as archive:
        members = (
            (
                member.name,
                lambda member=member: _require_tar_member(
                    archive,
                    member,
                ),
            )
            for member in archive.getmembers()
            if member.isfile()
        )
        _copy_release_members(
            members,
            expected_by_basename,
            directory,
            extracted,
            skipped_roles,
        )


def _copy_release_members(
    members: Iterator[tuple[str, Callable[[], BinaryIO]]],
    expected_by_basename: Mapping[str, str],
    directory: Path,
    extracted: dict[str, Path],
    skipped_roles: set[str],
) -> None:
    for member_name, opener in members:
        _validate_archive_member(member_name)
        base_name = PurePosixPath(member_name).name.lower()
        table_name = expected_by_basename.get(base_name)
        if table_name is None or table_name in skipped_roles:
            continue
        if table_name in extracted:
            raise ValueError(
                f"ChEBI release archive contains multiple candidates for "
                f"{table_name}: {member_name}"
            )
        suffix = ".tsv.gz" if base_name.endswith(".gz") else ".tsv"
        path = directory / f"{table_name}{suffix}"
        with opener() as handle_in, path.open("wb") as handle_out:
            shutil.copyfileobj(handle_in, handle_out)
        extracted[table_name] = path


def _require_tar_member(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
) -> BinaryIO:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"Cannot read archive member: {member.name}")
    return cast(BinaryIO, handle)


@contextmanager
def _open_ontology_text(path: Path) -> Generator[TextIO]:
    if _is_gzip(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            member = _select_ontology_member(
                member.filename for member in archive.infolist() if not member.is_dir()
            )
            with (
                archive.open(member) as raw,
                io.TextIOWrapper(
                    raw,
                    encoding="utf-8",
                    errors="replace",
                ) as handle,
            ):
                yield handle
        return
    if tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            members = {
                member.name: member
                for member in archive.getmembers()
                if member.isfile()
            }
            member_name = _select_ontology_member(iter(members))
            raw = archive.extractfile(members[member_name])
            if raw is None:
                raise ValueError(f"Cannot read ontology archive member: {member_name}")
            with (
                raw,
                io.TextIOWrapper(
                    raw,
                    encoding="utf-8",
                    errors="replace",
                ) as handle,
            ):
                yield handle
        return
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        yield handle


def _inspect_ontology_source(path: Path) -> None:
    with _open_ontology_text(path) as handle:
        for line in handle:
            if line.strip() == "[Term]":
                return
    raise ValueError(f"No [Term] stanza found in ontology source: {path}")


def _select_ontology_member(member_names: Iterator[str]) -> str:
    members = list(member_names)
    for member in members:
        _validate_archive_member(member)
    candidates = [
        member
        for member in members
        if PurePosixPath(member).name.lower().endswith(".obo")
    ]
    if len(candidates) != 1:
        raise ValueError(
            "Ontology archive must contain exactly one .obo member; "
            f"found {len(candidates)}"
        )
    return candidates[0]


def _validate_archive_member(member_name: str) -> None:
    path = PurePosixPath(member_name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe archive member path: {member_name}")


def _require_compounds(table_sources: Mapping[str, Path]) -> None:
    if "compound" not in table_sources:
        raise ValueError("ChEBI release does not contain compounds.tsv")


def _validate_table_source_inventory(table_sources: Mapping[str, Path]) -> None:
    physical_roles: dict[tuple[int, int], tuple[str, Path]] = {}
    for table_name, path in table_sources.items():
        resolved = path.resolve()
        stat = resolved.stat()
        physical_id = (stat.st_dev, stat.st_ino)
        if physical_id in physical_roles:
            other_name, other_path = physical_roles[physical_id]
            raise ValueError(
                f"ChEBI input roles {other_name!r} and {table_name!r} refer to "
                f"the same physical file: {other_path}"
            )
        physical_roles[physical_id] = (table_name, resolved)


def _require_file_or_directory(path: os.PathLike[str] | str) -> Path:
    candidate = Path(path)
    if not candidate.exists():
        raise FileNotFoundError(candidate)
    return candidate


def _require_file(path: os.PathLike[str] | str) -> Path:
    candidate = _require_file_or_directory(path)
    if not candidate.is_file():
        raise ValueError(f"Expected a file: {candidate}")
    return candidate


def _optional_file(path: os.PathLike[str] | str | None) -> Path | None:
    return None if path is None else _require_file(path)


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _source_record(logical_name: str, path: Path) -> SourceFileRecord:
    return SourceFileRecord(
        logical_name=logical_name,
        path=path,
        media_type=_media_type(path),
    )


def _media_type(path: Path) -> str:
    if _is_gzip(path):
        return "application/gzip"
    if zipfile.is_zipfile(path):
        return "application/zip"
    if tarfile.is_tarfile(path):
        return "application/x-tar"
    if path.name.lower().endswith(".obo"):
        return "text/obo"
    return "text/tab-separated-values"


def _scope(snapshot: _ChEBISnapshot, has_tables: bool) -> str:
    components: list[str] = []
    if has_tables:
        components.append("tables")
    if snapshot.file_obo is not None:
        components.append("chebi_ontology")
    if snapshot.file_chemont_obo is not None:
        components.append("chemont")
    return "+".join(components)


def _first(stanza: Mapping[str, list[str]], key: str) -> str | None:
    values = stanza.get(key)
    return values[0] if values else None


def _quoted(value: str | None) -> str | None:
    if value is None or not value.startswith('"'):
        return value
    escaped = False
    result: list[str] = []
    for character in value[1:]:
        if character == '"' and not escaped:
            return "".join(result)
        if character == "\\" and not escaped:
            escaped = True
            continue
        result.append(character)
        escaped = False
    return "".join(result)


def _synonym_scope(value: str) -> str:
    remainder = value
    if value.startswith('"'):
        escaped = False
        for index, character in enumerate(value[1:], start=1):
            if character == '"' and not escaped:
                remainder = value[index + 1 :].strip()
                break
            escaped = character == "\\" and not escaped
    return remainder.split(maxsplit=1)[0] if remainder else ""
