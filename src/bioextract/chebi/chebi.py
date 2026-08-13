from __future__ import annotations

import gzip
import os
import shutil
import stat
import tarfile
import tempfile
import zipfile
from collections.abc import Generator, Iterable, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Literal, TextIO

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
    archive_table_members: Mapping[str, str] = field(default_factory=dict[str, str])
    archive_obo_member: str | None = None
    archive_sdf_member: str | None = None
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

        ``source`` may be a table-release directory, zip archive, or tar
        archive. Explicit non-``None`` roles replace discovered roles; without
        a source, ``compounds`` remains the required primary role. OBO files
        are never selected by this representation-specific constructor.

        Examples:
            Replace one discovered table role while retaining the source
            profile:

            >>> db = ChEBIDatabase.from_table_files(  # doctest: +SKIP
            ...     "chebi-release",
            ...     names="overrides/names.tsv",
            ... )

            Missing the required primary table fails before publication:

            >>> try:
            ...     ChEBIDatabase.from_table_files()
            ... except ValueError as error:
            ...     print(error)
            ChEBI release does not contain compounds.tsv
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
        chemont = _optional_file(chemont_obo)
        _validate_physical_inventory(
            {
                **explicit_sources,
                **({"chemont_obo": chemont} if chemont is not None else {}),
            }
        )
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
                skipped_roles=set(explicit_sources),
            )
            table_sources = {**discovered, **explicit_sources}
            _require_compounds(table_sources)
            _validate_physical_inventory(
                {
                    **table_sources,
                    **({"chemont_obo": chemont} if chemont is not None else {}),
                }
            )
            return cls(
                snapshot=_ChEBISnapshot(
                    table_sources=table_sources,
                    file_chemont_obo=chemont,
                )
            )

        archive_members = _inspect_table_archive(
            source_path,
            skipped_roles=set(explicit_sources),
        )
        if "compound" not in archive_members and "compound" not in explicit_sources:
            raise ValueError("ChEBI release does not contain compounds.tsv")
        _validate_physical_inventory(
            {
                "release_archive": source_path,
                **explicit_sources,
                **({"chemont_obo": chemont} if chemont is not None else {}),
            }
        )
        return cls(
            snapshot=_ChEBISnapshot(
                table_sources=explicit_sources,
                release_source=source_path,
                archive_table_members=archive_members,
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

        ``source`` may be an exact OBO file, a directory containing exactly one
        OBO candidate, or an archive containing exactly one OBO member. Only
        directory and archive sources discover a matching SDF supplement; an
        exact file requires an explicit ``sdf``.

        Examples:
            A missing ontology input fails at construction:

            >>> try:
            ...     ChEBIDatabase.from_obo("missing-chebi.obo")
            ... except FileNotFoundError as error:
            ...     print(error.filename)
            missing-chebi.obo
        """
        chemont = _optional_file(chemont_obo)
        source_path = _require_file_or_directory(source)
        explicit_sdf = _optional_file(sdf)
        if source_path.is_dir():
            file_obo = _discover_ontology_file(
                source_path,
                excluded=(() if chemont is None else (chemont,)),
            )
            file_sdf = (
                explicit_sdf
                if explicit_sdf is not None
                else _discover_supplement_file(source_path, suffix=".sdf")
            )
            _validate_physical_inventory(
                {
                    "chebi_obo": file_obo,
                    **({"chebi_sdf": file_sdf} if file_sdf is not None else {}),
                    **({"chemont_obo": chemont} if chemont is not None else {}),
                }
            )
            _inspect_ontology_source(file_obo)
            if file_sdf is not None:
                _inspect_sdf_source(file_sdf)
            return cls(
                snapshot=_ChEBISnapshot(
                    file_obo=file_obo,
                    file_sdf=file_sdf,
                    file_chemont_obo=chemont,
                )
            )

        if _is_archive(source_path):
            archive_members = _inspect_ontology_archive(
                source_path,
                discover_sdf=explicit_sdf is None,
            )
            obo_member = archive_members["obo"]
            sdf_member = archive_members.get("sdf")
            _validate_physical_inventory(
                {
                    "release_archive": source_path,
                    **({"chebi_sdf": explicit_sdf} if explicit_sdf is not None else {}),
                    **({"chemont_obo": chemont} if chemont is not None else {}),
                }
            )
            _inspect_archive_member(source_path, obo_member, suffix=".obo")
            if sdf_member is not None:
                _inspect_archive_member(source_path, sdf_member, suffix=".sdf")
            if explicit_sdf is not None:
                _inspect_sdf_source(explicit_sdf)
            if chemont is not None:
                _inspect_ontology_source(chemont)
            return cls(
                snapshot=_ChEBISnapshot(
                    release_source=source_path,
                    archive_obo_member=obo_member,
                    archive_sdf_member=(
                        None if explicit_sdf is not None else sdf_member
                    ),
                    file_sdf=explicit_sdf,
                    file_chemont_obo=chemont,
                )
            )

        _validate_physical_inventory(
            {
                "chebi_obo": source_path,
                **({"chebi_sdf": explicit_sdf} if explicit_sdf is not None else {}),
                **({"chemont_obo": chemont} if chemont is not None else {}),
            }
        )
        _inspect_ontology_source(source_path)
        if explicit_sdf is not None:
            _inspect_sdf_source(explicit_sdf)
        if chemont is not None:
            _inspect_ontology_source(chemont)
        return cls(
            snapshot=_ChEBISnapshot(
                file_obo=source_path,
                file_sdf=explicit_sdf,
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
            >>> selection.matches().select("chebi_id").collect().unique().to_list()  # doctest: +SKIP
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
            >>> selection.matches().select("group_id").collect().to_list()  # doctest: +SKIP
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

        with _prepared_sources(self.snapshot) as (
            table_sources,
            file_obo,
            file_sdf,
            file_chemont_obo,
            provenance_sources,
        ):
            relations: list[RelationSpec] = []
            validation_issues = ()
            if file_obo is not None:
                canonical = build_canonical_relations(
                    file_obo,
                    file_sdf=file_sdf,
                )
                relations.extend(canonical.relations)
                validation_issues = canonical.validation_issues
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
            if file_chemont_obo is not None:
                relations.extend(
                    _obo_relations(
                        file_chemont_obo,
                        namespace="chemont",
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
                scope=_scope(
                    has_tables=bool(table_sources),
                    has_obo=file_obo is not None,
                    has_chemont=file_chemont_obo is not None,
                ),
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
def _prepared_sources(
    snapshot: _ChEBISnapshot,
) -> Generator[
    tuple[
        dict[str, Path], Path | None, Path | None, Path | None, list[SourceFileRecord]
    ]
]:
    """Materialize archive members and normalize compressed table inputs."""
    with tempfile.TemporaryDirectory(prefix="bioextract-chebi-") as value:
        directory = Path(value)
        table_sources = dict(snapshot.table_sources)
        file_obo = snapshot.file_obo
        file_sdf = snapshot.file_sdf
        provenance_sdf = file_sdf
        file_chemont_obo = snapshot.file_chemont_obo
        provenance: list[SourceFileRecord] = []

        archive = snapshot.release_source
        archive_used = False
        if archive is not None:
            for table_name, member_name in snapshot.archive_table_members.items():
                table_sources[table_name] = _extract_archive_member(
                    archive,
                    member_name,
                    directory / f"{table_name}.source",
                )
                archive_used = True
            if snapshot.archive_obo_member is not None:
                file_obo = _extract_archive_member(
                    archive,
                    snapshot.archive_obo_member,
                    directory / "chebi.obo.source",
                )
                archive_used = True
            if snapshot.archive_sdf_member is not None and file_sdf is None:
                file_sdf = _extract_archive_member(
                    archive,
                    snapshot.archive_sdf_member,
                    directory / "chebi.sdf.source",
                )
                archive_used = True
            if archive_used:
                provenance.append(_source_record("release_archive", archive))

        prepared_tables: dict[str, Path] = {}
        for table_name, file_source in table_sources.items():
            if _is_gzip(file_source):
                file_plain = directory / f"{table_name}.tsv"
                with (
                    gzip.open(file_source, "rb") as handle_in,
                    file_plain.open("wb") as handle_out,
                ):
                    shutil.copyfileobj(handle_in, handle_out)
                prepared_tables[table_name] = file_plain
            else:
                prepared_tables[table_name] = file_source

        if file_sdf is not None and _is_archive(file_sdf):
            sdf_member = _one_candidate(
                "Explicit SDF archive must contain exactly one .sdf member",
                [
                    member
                    for member in _archive_member_names(file_sdf)
                    if _normalized_suffix_name(
                        PurePosixPath(member).name,
                        ".sdf",
                    )
                ],
            )
            file_sdf = _extract_archive_member(
                file_sdf,
                sdf_member,
                directory / "explicit-chebi.sdf.source",
            )

        provenance.extend(
            _source_record(table_name, file_source)
            for table_name, file_source in snapshot.table_sources.items()
        )
        if file_obo is not None and snapshot.archive_obo_member is None:
            provenance.append(_source_record("chebi_obo", file_obo))
        if provenance_sdf is not None and snapshot.archive_sdf_member is None:
            provenance.append(_source_record("chebi_sdf", provenance_sdf))
        if file_chemont_obo is not None:
            provenance.append(_source_record("chemont_obo", file_chemont_obo))
        yield (
            prepared_tables,
            file_obo,
            file_sdf,
            file_chemont_obo,
            provenance,
        )


def _discover_table_files(
    directory: Path,
    *,
    skipped_roles: set[str] | None = None,
) -> dict[str, Path]:
    skipped: set[str] = set() if skipped_roles is None else skipped_roles
    candidates: dict[str, list[Path]] = {table_name: [] for table_name in _TABLE_FILES}
    for candidate in sorted(directory.rglob("*")):
        if not candidate.is_file():
            continue
        table_name = _table_role_for_name(candidate.name)
        if table_name is not None and table_name not in skipped:
            candidates[table_name].append(candidate)
    return {
        table_name: _one_candidate(
            f"ChEBI table role {table_name} must contain exactly one candidate",
            paths,
        )
        for table_name, paths in candidates.items()
        if paths
    }


def _discover_ontology_file(
    directory: Path,
    *,
    excluded: Iterable[Path] = (),
) -> Path:
    excluded_paths = tuple(excluded)
    candidates = [
        candidate
        for candidate in sorted(directory.rglob("*"))
        if candidate.is_file()
        and _normalized_suffix_name(candidate.name, ".obo")
        and not any(_same_physical(candidate, path) for path in excluded_paths)
    ]
    return _one_candidate(
        "ChEBI source directory must contain exactly one .obo candidate",
        candidates,
    )


def _discover_supplement_file(directory: Path, *, suffix: str) -> Path | None:
    candidates = [
        candidate
        for candidate in sorted(directory.rglob("*"))
        if candidate.is_file() and _normalized_suffix_name(candidate.name, suffix)
    ]
    return (
        None
        if not candidates
        else _one_candidate(
            f"ChEBI source directory must contain exactly one {suffix} candidate",
            candidates,
        )
    )


def _inspect_table_archive(
    path: Path,
    *,
    skipped_roles: set[str] | None = None,
) -> dict[str, str]:
    skipped: set[str] = set() if skipped_roles is None else skipped_roles
    members = _archive_member_names(path)
    candidates: dict[str, list[str]] = {table_name: [] for table_name in _TABLE_FILES}
    for member_name in members:
        table_name = _table_role_for_name(PurePosixPath(member_name).name)
        if table_name is not None and table_name not in skipped:
            candidates[table_name].append(member_name)
    return {
        table_name: _one_candidate(
            f"ChEBI archive table role {table_name} must contain exactly one candidate",
            names,
        )
        for table_name, names in candidates.items()
        if names
    }


def _inspect_ontology_archive(
    path: Path,
    *,
    discover_sdf: bool = True,
) -> dict[str, str]:
    members = _archive_member_names(path)
    obo = _one_candidate(
        "ChEBI source archive must contain exactly one .obo member",
        [
            member
            for member in members
            if _normalized_suffix_name(PurePosixPath(member).name, ".obo")
        ],
    )
    sdf_candidates = (
        [
            member
            for member in members
            if _normalized_suffix_name(PurePosixPath(member).name, ".sdf")
        ]
        if discover_sdf
        else []
    )
    return {
        "obo": obo,
        **(
            {
                "sdf": _one_candidate(
                    "ChEBI source archive must contain exactly one .sdf member",
                    sdf_candidates,
                )
            }
            if sdf_candidates
            else {}
        ),
    }


def _inspect_archive_member(path: Path, member_name: str, *, suffix: str) -> None:
    with tempfile.TemporaryDirectory(prefix="bioextract-chebi-inspect-") as value:
        extracted = _extract_archive_member(path, member_name, Path(value) / "member")
        if suffix == ".obo":
            _inspect_ontology_source(extracted)
        elif suffix == ".sdf":
            _inspect_sdf_source(extracted)
        else:
            raise ValueError(f"Unsupported archive inspection suffix: {suffix}")


def _archive_member_names(path: Path) -> tuple[str, ...]:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            names: list[str] = []
            for info in archive.infolist():
                _validate_archive_member(info.filename)
                if info.is_dir():
                    continue
                mode = (info.external_attr >> 16) & 0xFFFF
                if info.create_system == 3 and mode and stat.S_ISLNK(mode):
                    raise ValueError(f"Archive symlink is not allowed: {info.filename}")
                names.append(_normalize_member_name(info.filename))
    elif tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            names = []
            for member in archive.getmembers():
                _validate_archive_member(member.name)
                if member.issym() or member.islnk():
                    raise ValueError(f"Archive link is not allowed: {member.name}")
                if member.isfile():
                    names.append(_normalize_member_name(member.name))
    else:
        raise ValueError(
            f"ChEBI source is not a directory, zip, or tar archive: {path}"
        )
    if len(set(names)) != len(names):
        raise ValueError(f"Archive contains duplicate member paths: {path}")
    return tuple(names)


def _extract_archive_member(path: Path, member_name: str, destination: Path) -> Path:
    normalized = _normalize_member_name(member_name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            matches = [
                info
                for info in archive.infolist()
                if not info.is_dir()
                and _normalize_member_name(info.filename) == normalized
            ]
            if len(matches) != 1:
                raise ValueError(f"Archive member is not unique: {member_name}")
            with (
                archive.open(matches[0]) as handle_in,
                destination.open("wb") as handle_out,
            ):
                shutil.copyfileobj(handle_in, handle_out)
    elif tarfile.is_tarfile(path):
        with tarfile.open(path, mode="r:*") as archive:
            matches = [
                member
                for member in archive.getmembers()
                if member.isfile() and _normalize_member_name(member.name) == normalized
            ]
            if len(matches) != 1:
                raise ValueError(f"Archive member is not unique: {member_name}")
            handle_in = archive.extractfile(matches[0])
            if handle_in is None:
                raise ValueError(f"Cannot read archive member: {member_name}")
            with handle_in, destination.open("wb") as handle_out:
                shutil.copyfileobj(handle_in, handle_out)
    else:
        raise ValueError(f"Not an archive: {path}")
    return destination


@contextmanager
def _open_ontology_text(path: Path) -> Generator[TextIO]:
    with _open_text_source(path, preferred_suffix=".obo") as handle:
        yield handle


@contextmanager
def _open_text_source(path: Path, *, preferred_suffix: str) -> Generator[TextIO]:
    if _is_gzip(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return
    if _is_archive(path):
        members = _archive_member_names(path)
        member_name = _one_candidate(
            f"Archive must contain exactly one {preferred_suffix} member",
            [
                member
                for member in members
                if _normalized_suffix_name(
                    PurePosixPath(member).name,
                    preferred_suffix,
                )
            ],
        )
        with tempfile.TemporaryDirectory(prefix="bioextract-chebi-text-") as value:
            extracted = _extract_archive_member(
                path, member_name, Path(value) / "source"
            )
            with _open_text_source(
                extracted, preferred_suffix=preferred_suffix
            ) as handle:
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


def _inspect_sdf_source(path: Path) -> None:
    with _open_text_source(path, preferred_suffix=".sdf") as handle:
        if any(line.strip() == "$$$$" for line in handle):
            return
    raise ValueError(f"No SDF record found in structure source: {path}")


def _table_role_for_name(name: str) -> str | None:
    normalized = name.lower()
    if normalized.endswith(".gz"):
        normalized = normalized[:-3]
    for table_name, base_name in _TABLE_FILES.items():
        if normalized == base_name.lower():
            return table_name
    return None


def _normalized_suffix_name(name: str, suffix: str) -> bool:
    normalized = name.lower()
    if normalized.endswith(".gz"):
        normalized = normalized[:-3]
    return normalized.endswith(suffix.lower())


def _one_candidate[T](label: str, candidates: Iterable[T]) -> T:
    values = list(candidates)
    if len(values) != 1:
        detail = ", ".join(str(value) for value in values)
        raise ValueError(f"{label}; found {len(values)}: {detail}")
    return values[0]


def _normalize_member_name(member_name: str) -> str:
    _validate_archive_member(member_name)
    return PurePosixPath(member_name).as_posix()


def _validate_archive_member(member_name: str) -> None:
    path = PurePosixPath(member_name)
    if not member_name or path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe archive member path: {member_name}")


def _validate_physical_inventory(sources: Mapping[str, Path]) -> None:
    identities: dict[tuple[int, int], str] = {}
    for logical_name, path in sources.items():
        identity = _physical_identity(path)
        previous = identities.get(identity)
        if previous is not None:
            raise ValueError(
                "ChEBI source roles must use different physical files: "
                f"{previous} and {logical_name} both refer to {path}"
            )
        identities[identity] = logical_name


def _physical_identity(path: Path) -> tuple[int, int]:
    stat_result = path.stat()
    return stat_result.st_dev, stat_result.st_ino


def _same_physical(left: Path, right: Path) -> bool:
    try:
        return _physical_identity(left) == _physical_identity(right)
    except FileNotFoundError:
        return False


def _require_compounds(table_sources: Mapping[str, Path]) -> None:
    if "compound" not in table_sources:
        raise ValueError("ChEBI release does not contain compounds.tsv")


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


def _is_archive(path: Path) -> bool:
    return zipfile.is_zipfile(path) or tarfile.is_tarfile(path)


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


def _scope(*, has_tables: bool, has_obo: bool, has_chemont: bool) -> str:
    components: list[str] = []
    if has_tables:
        components.append("tables")
    if has_obo:
        components.append("chebi_ontology")
    if has_chemont:
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
