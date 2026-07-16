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
        Disable the interactions-file size guard while keeping other defaults:

        >>> limits = OmniPathResourceLimits(file_interactions_bytes_max=None)
        >>> limits.file_interactions_bytes_max is None
        True
    """

    file_enzsub_bytes_max: int | None = 512 * 1024 * 1024
    file_interactions_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000
