from __future__ import annotations

import gzip
import inspect
import warnings
from pathlib import Path

from bioextract.eggnog import EggNOGDatabase
from bioextract.eggnog.util import is_gzip_file

WARNING_MESSAGE = (
    "Compressed eggNOG SQLite source detected. bioextract must fully "
    "decompress it to a temporary SQLite file before access, and the "
    "decompressed file is not persisted. Repeated use may repeat this cost; "
    "for long-term use, decompress the source once and pass the .db file."
)


def test_gzip_transport_detection_uses_magic_bytes(tmp_path: Path) -> None:
    wrapped_without_suffix = tmp_path / "eggnog.db"
    with gzip.open(wrapped_without_suffix, "wb") as handle:
        handle.write(b"SQLite bytes")
    plain_with_suffix = tmp_path / "plain.gz"
    plain_with_suffix.write_bytes(b"SQLite format 3\x00")

    assert is_gzip_file(wrapped_without_suffix)
    assert not is_gzip_file(plain_with_suffix)


def test_gzip_warning_has_exact_message_and_caller_location(tmp_path: Path) -> None:
    source = tmp_path / "eggnog.db.gz"
    with gzip.open(source, "wb") as handle:
        handle.write(b"SQLite bytes")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        caller_line = inspect.currentframe().f_lineno + 1  # type: ignore[union-attr]
        EggNOGDatabase.from_sqlite(source)

    assert len(caught) == 1
    assert caught[0].category is UserWarning
    assert str(caught[0].message) == WARNING_MESSAGE
    assert caught[0].filename == __file__
    assert caught[0].lineno == caller_line
