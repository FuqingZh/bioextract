import os

import duckdb
import polars as pl
import pytest

from bioextract._runtime import (
    DEFAULT_TEST_THREADS,
    NATIVE_THREAD_ENV_VARS,
    configure_validation_environment,
)


def test_unset_validation_budget_uses_context_default() -> None:
    environment: dict[str, str] = {}

    actual = configure_validation_environment(environ=environment)

    assert actual == DEFAULT_TEST_THREADS
    assert environment["BIOEXTRACT_TEST_THREADS"] == str(DEFAULT_TEST_THREADS)


@pytest.mark.parametrize("budget", ["1", "4", "127"])
def test_validation_budget_accepts_every_positive_integer(budget: str) -> None:
    environment = {"BIOEXTRACT_TEST_THREADS": budget}

    actual = configure_validation_environment(environ=environment)

    assert actual == int(budget)
    assert environment["BIOEXTRACT_TEST_THREADS"] == budget


@pytest.mark.parametrize("invalid", ["0", "-1", "1.5", "many"])
def test_validation_budget_rejects_non_positive_integer(invalid: str) -> None:
    with pytest.raises(ValueError, match="must be a positive integer"):
        configure_validation_environment(environ={"BIOEXTRACT_TEST_THREADS": invalid})


@pytest.mark.parametrize("variable", NATIVE_THREAD_ENV_VARS)
@pytest.mark.parametrize(
    ("caller_value", "expected"),
    [(None, "4"), ("2", "2"), ("4", "4"), ("12", "4")],
)
def test_each_native_pool_is_set_preserved_or_clamped(
    variable: str, caller_value: str | None, expected: str
) -> None:
    environment = {"BIOEXTRACT_TEST_THREADS": "4"}
    if caller_value is not None:
        environment[variable] = caller_value

    configure_validation_environment(environ=environment)

    assert environment[variable] == expected


@pytest.mark.parametrize("variable", NATIVE_THREAD_ENV_VARS)
@pytest.mark.parametrize("invalid", ["0", "-1", "1.5", "many"])
def test_each_native_pool_rejects_invalid_value(variable: str, invalid: str) -> None:
    with pytest.raises(ValueError, match=variable):
        configure_validation_environment(
            environ={"BIOEXTRACT_TEST_THREADS": "4", variable: invalid}
        )


def test_analytical_engines_use_bounded_test_threads() -> None:
    expected = int(os.environ["BIOEXTRACT_TEST_THREADS"])

    assert pl.thread_pool_size() <= expected
    default_threads = duckdb.sql(
        "SELECT current_setting('threads')::INTEGER"
    ).fetchone()
    assert default_threads == (expected,)
    with duckdb.connect() as connection:
        actual = connection.execute(
            "SELECT current_setting('threads')::INTEGER"
        ).fetchone()
    assert actual == (expected,)
