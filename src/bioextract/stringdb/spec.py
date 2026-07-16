from dataclasses import dataclass

__all__ = [
    "StringResourceLimits",
]


@dataclass(frozen=True, slots=True)
class StringResourceLimits:
    """Fail-fast limits for STRING snapshot files and normalized selections.

    File limits are measured in bytes on disk. Input and group limits are
    applied after whitespace cleanup, UniProt pipe-ID normalization, blank
    removal, and deduplication. Set any field to ``None`` to disable that
    check.

    Examples:
        Reject an oversized alias snapshot before parsing it:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> from bioextract.stringdb import StringDb
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_aliases = Path(dir_tmp) / "9606.protein.aliases.v12.0.txt"
        ...     _ = file_aliases.write_text(
        ...         "string_protein_id\\talias\\tsource\\n"
        ...     )
        ...     limits = StringResourceLimits(file_aliases_bytes_max=1)
        ...     try:
        ...         StringDb.from_files(file_aliases=file_aliases, limits=limits)
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_aliases_bytes_max: int | None = 512 * 1024 * 1024
    file_links_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000
