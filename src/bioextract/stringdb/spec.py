from dataclasses import dataclass

__all__ = [
    "StringResourceLimits",
]


@dataclass(frozen=True, slots=True)
class StringResourceLimits:
    file_aliases_bytes_max: int | None = 512 * 1024 * 1024
    file_links_bytes_max: int | None = 4 * 1024 * 1024 * 1024
    num_input_ids_max: int | None = 100_000
    num_groups_max: int | None = 1_000
