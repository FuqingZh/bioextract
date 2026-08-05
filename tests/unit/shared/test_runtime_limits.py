import os
from typing import Any, cast

import duckdb
import numpy as np
import polars as pl
from threadpoolctl import threadpool_info  # pyright: ignore[reportMissingTypeStubs]


def test_analytical_engines_use_bounded_test_threads() -> None:
    expected = int(os.environ["BIOEXTRACT_TEST_THREADS"])

    assert pl.thread_pool_size() <= expected
    with duckdb.connect() as connection:
        actual = connection.execute(
            "SELECT current_setting('threads')::INTEGER"
        ).fetchone()
    assert actual == (expected,)

    np.dot(np.ones((32, 32)), np.ones((32, 32)))
    for pool in cast(list[dict[str, Any]], threadpool_info()):
        thread_count = pool.get("num_threads")
        if isinstance(thread_count, int):
            assert thread_count <= expected, pool
