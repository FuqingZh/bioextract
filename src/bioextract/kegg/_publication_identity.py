from __future__ import annotations

from pathlib import Path

from bioextract.errors import IntegrityError

type FileIdentity = tuple[int, int, int, int, int]


def capture_file_identity(path: Path) -> FileIdentity:
    """Capture the filesystem identity used to pin one KEGG publication."""
    current = path.stat()
    return (
        current.st_dev,
        current.st_ino,
        current.st_size,
        current.st_mtime_ns,
        current.st_ctime_ns,
    )


def assert_publication_current(
    path: Path,
    expected: FileIdentity | None,
    *,
    profile_label: str,
) -> None:
    """Reject a reopened path that no longer names the validated publication."""
    if expected is None:
        return
    try:
        current = capture_file_identity(path)
    except OSError as error:
        raise IntegrityError(
            f"KEGG {profile_label} publication was replaced; "
            "reopen it with from_duckdb()"
        ) from error
    if current != expected:
        raise IntegrityError(
            f"KEGG {profile_label} publication was replaced; "
            "reopen it with from_duckdb()"
        )
