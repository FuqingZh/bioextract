"""Direct, validated access to official curated ComplexTab snapshots."""

from dataclasses import dataclass, replace
from pathlib import Path

import polars as pl

__all__ = ["ComplexPortalDatabase"]

COLUMNS = (
    "#Complex ac",
    "Recommended name",
    "Aliases for complex",
    "Taxonomy identifier",
    "Identifiers (and stoichiometry) of molecules in complex",
    "Evidence Code",
    "Experimental evidence",
    "Go Annotations",
    "Cross references",
    "Description",
    "Complex properties",
    "Complex assembly",
    "Ligand",
    "Disease",
    "Agonist",
    "Antagonist",
    "Comment",
    "Source",
    "Expanded participant list",
)
_ACCESSION = r"[A-Z0-9]{6,10}(?:-(?:[0-9]+|PRO_[0-9]+))?"
_MEMBER = rf"^({_ACCESSION}|EBI-[0-9]+|\[{_ACCESSION}(?:,{_ACCESSION})*\])\([0-9]+\)$"


def _identity(path: Path) -> tuple[int, int, int, int]:
    info = path.stat()
    return info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def _scan(path: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        path,
        separator="\t",
        schema=dict.fromkeys(COLUMNS, pl.String),
        null_values=["-", ""],
        quote_char=None,
    )


@dataclass(frozen=True, slots=True)
class ComplexPortalDatabase:
    """A caller-owned official ComplexTab file with explicit release identity.

    Use ``from_files`` to validate the snapshot. Native lazy relations are
    replayable; no download, gene mapping or analysis policy is performed here.
    """

    path: Path
    release_version: str
    _identity: tuple[int, int, int, int]
    taxon_id: str | None = None

    @classmethod
    def from_files(
        cls, *, complextab: Path | str, release_version: str
    ) -> "ComplexPortalDatabase":
        """Validate official headers, complex identities and expanded member tokens.

        A missing or empty source is an error. The caller must supply the curated
        product and retain acquisition provenance; filenames never set identity.
        """
        path = Path(complextab).resolve(strict=True)
        if not release_version.strip() or release_version != release_version.strip():
            raise ValueError("Complex Portal release_version must be explicit")
        before = _identity(path)
        with path.open(encoding="utf-8") as source:
            if tuple(source.readline().rstrip("\r\n").split("\t")) != COLUMNS:
                raise ValueError(
                    "ComplexTab requires the official ordered 19-column header"
                )
        table = _scan(path)
        identities = table.select(
            "#Complex ac",
            "Recommended name",
            "Taxonomy identifier",
            "Evidence Code",
            "Expanded participant list",
        ).collect()
        if identities.is_empty() or any(
            identities[name].null_count() for name in identities.columns
        ):
            raise ValueError(
                "ComplexTab contains missing required identities or members"
            )
        if identities["#Complex ac"].is_duplicated().any():
            raise ValueError("duplicate Complex Portal accession")
        if identities.filter(
            ~pl.col("#Complex ac").str.contains(r"^CPX-[0-9]+$")
            | ~pl.col("Taxonomy identifier").str.contains(r"^[1-9][0-9]*$")
        ).height:
            raise ValueError("invalid Complex Portal accession or taxon")
        tokens = identities.select(
            "#Complex ac",
            pl.col("Expanded participant list").str.split("|").alias("token"),
        ).explode("token")
        if tokens.filter(~pl.col("token").str.contains(_MEMBER)).height:
            raise ValueError(
                "invalid expanded ComplexTab protein/set/stoichiometry token"
            )
        if tokens.select("#Complex ac", "token").is_duplicated().any():
            raise ValueError("duplicate expanded protein member token")
        if _identity(path) != before:
            raise ValueError("ComplexTab mutated during validation")
        return cls(path, release_version, before)

    def _validated_scan(self) -> pl.LazyFrame:
        if _identity(self.path) != self._identity:
            raise ValueError("ComplexTab changed; reopen the source")
        frame = _scan(self.path)
        return (
            frame
            if self.taxon_id is None
            else frame.filter(pl.col("Taxonomy identifier") == self.taxon_id)
        )

    def with_taxon(self, taxon_id: str) -> "ComplexPortalDatabase":
        """Return a taxon-scoped view; no species-name or orthology inference."""
        if not taxon_id.isdecimal() or int(taxon_id) < 1:
            raise ValueError("Complex Portal taxon must be a positive NCBI identifier")
        return replace(self, taxon_id=taxon_id)

    def complexes(self) -> pl.LazyFrame:
        """Return the unchanged official complex-level columns in source order."""
        return self._validated_scan()

    def expanded_members(self) -> pl.LazyFrame:
        """Expand official participant identities and alternative sets with exact suffixes.

        Alternatives share ``member_group_id``; they do not assert simultaneous
        membership. Unknown stoichiometry (0) becomes null while the raw token
        remains available. Nested complexes use the official expanded list.
        """
        frame = (
            self._validated_scan()
            .select(
                pl.col("#Complex ac").alias("complex_id"),
                pl.col("Recommended name").alias("complex_name"),
                pl.col("Taxonomy identifier").alias("taxon_id"),
                pl.col("Evidence Code").alias("evidence_code"),
                pl.col("Expanded participant list")
                .str.split("|")
                .alias("member_token"),
            )
            .explode("member_token")
            .with_columns(
                pl.col("member_token")
                .str.extract(r"^(.+)\([0-9]+\)$", 1)
                .alias("member_expression"),
                pl.col("member_token")
                .str.extract(r"\(([0-9]+)\)$", 1)
                .cast(pl.Int64)
                .alias("stoichiometry"),
            )
            .with_columns(
                pl.concat_str("complex_id", "member_token", separator=":").alias(
                    "member_group_id"
                ),
                pl.col("member_expression")
                .str.starts_with("[")
                .alias("is_alternative"),
                pl.col("member_expression")
                .str.strip_chars("[]")
                .str.split(",")
                .alias("participant_id"),
                pl.when(pl.col("stoichiometry") == 0)
                .then(None)
                .otherwise(pl.col("stoichiometry"))
                .alias("stoichiometry"),
                pl.lit(self.release_version).alias("resource_version"),
            )
            .explode("participant_id")
        )
        frame = frame.with_columns(
            pl.when(pl.col("participant_id").str.starts_with("EBI-"))
            .then(pl.lit("intact"))
            .otherwise(pl.lit("uniprot"))
            .alias("identifier_namespace")
        )
        return frame.select(
            "complex_id",
            "complex_name",
            "taxon_id",
            "member_group_id",
            "participant_id",
            "identifier_namespace",
            "is_alternative",
            "stoichiometry",
            "member_token",
            "evidence_code",
            "resource_version",
        )
