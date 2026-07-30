from dataclasses import dataclass

__all__ = [
    "OmniPathResourceLimits",
]


@dataclass(frozen=True, slots=True)
class OmniPathResourceLimits:
    """Fail-fast limits for OmniPath files and normalized selections.

    File limits are measured in bytes on disk. Input and group limits are
    applied after identifier and group normalization. Set any field to
    ``None`` to disable that check.

    Examples:
        Reject an oversized interaction snapshot before parsing it:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> from bioextract.omnipath import OmniPathDb
        >>> with TemporaryDirectory() as dir_tmp:
        ...     file_interactions = Path(dir_tmp) / "interactions.tsv"
        ...     _ = file_interactions.write_text("source\\ttarget\\n")
        ...     limits = OmniPathResourceLimits(file_interactions_bytes_max=1)
        ...     try:
        ...         OmniPathDb.from_files(
        ...             file_interactions=file_interactions,
        ...             limits=limits,
        ...         )
        ...     except ValueError as error:
        ...         print("exceeds configured size limit" in str(error))
        True
    """

    file_enzsub_bytes_max: int | None = 512 * 1024 * 1024
    file_interactions_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000
