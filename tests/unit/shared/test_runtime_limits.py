import os

import duckdb
import polars as pl


def test_analytical_engines_use_bounded_test_threads() -> None:
    expected = int(os.environ["BIOEXTRACT_TEST_THREADS"])

    assert pl.thread_pool_size() <= expected
    with duckdb.connect() as connection:
        actual = connection.execute(
            "SELECT current_setting('threads')::INTEGER"
        ).fetchone()
    assert actual == (expected,)
