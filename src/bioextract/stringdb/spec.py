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
        Limit one STRING query to 500 normalized identifiers:

        >>> limits = StringResourceLimits(num_input_ids_max=500)
        >>> limits.num_input_ids_max
        500
    """

    file_aliases_bytes_max: int | None = 512 * 1024 * 1024
    file_links_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000
