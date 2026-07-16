from __future__ import annotations

import doctest
import inspect
import sys
from collections.abc import Iterator

import pytest

import bioextract.go.ontology as go_ontology
from bioextract._tidy import TidyDataset
from bioextract.eggnog import EggnogDb, EggnogResourceLimits, EggnogSelection
from bioextract.go import GoDb, GoResourceLimits, GoSubsetId
from bioextract.go.ontology import run_tidy_go_ontology
from bioextract.interpro import (
    InterProDb,
    InterProResourceLimits,
    InterProSelection,
)
from bioextract.kegg import KeggDb, KeggResourceLimits
from bioextract.kegg.brite import (
    build_tidy_frames as build_kegg_tidy_frames,
    run_tidy_kegg_brite,
)
from bioextract.kegg.kegg import KeggSelection
from bioextract.omnipath import OmniPathDb, OmniPathResourceLimits
from bioextract.omnipath.omnipath import OmniPathSelection
from bioextract.reactome import ReactomeDb, ReactomeResourceLimits
from bioextract.reactome.reactome import ReactomeSelection
from bioextract.stringdb import StringDb, StringResourceLimits
from bioextract.stringdb.stringdb import StringSelection
from bioextract.uniprot import UniprotDb, UniprotResourceLimits
from bioextract.wikipathways import WikiPathwaysDb, WikiPathwaysResourceLimits
from bioextract.wikipathways.wikipathways import WikiPathwaysSelection

PUBLIC_CLASSES = (
    GoDb,
    GoResourceLimits,
    GoSubsetId,
    KeggDb,
    KeggResourceLimits,
    KeggSelection,
    ReactomeDb,
    ReactomeResourceLimits,
    ReactomeSelection,
    WikiPathwaysDb,
    WikiPathwaysResourceLimits,
    WikiPathwaysSelection,
    EggnogDb,
    EggnogResourceLimits,
    EggnogSelection,
    InterProDb,
    InterProResourceLimits,
    InterProSelection,
    UniprotDb,
    UniprotResourceLimits,
    StringDb,
    StringResourceLimits,
    StringSelection,
    OmniPathDb,
    OmniPathResourceLimits,
    OmniPathSelection,
    TidyDataset,
)

PUBLIC_FUNCTIONS = (
    ("go.ontology.run_tidy_go_ontology", run_tidy_go_ontology),
    ("kegg.brite.build_tidy_frames", build_kegg_tidy_frames),
    ("kegg.brite.run_tidy_kegg_brite", run_tidy_kegg_brite),
)

EXECUTABLE_DOCSTRING_TARGETS = (
    GoResourceLimits,
    KeggResourceLimits,
    ReactomeResourceLimits,
    WikiPathwaysResourceLimits,
    EggnogResourceLimits,
    InterProResourceLimits,
    UniprotResourceLimits,
    StringResourceLimits,
    OmniPathResourceLimits,
)

EXPECTED_PUBLIC_TARGET_COUNT = 124


def iter_public_docstring_targets() -> Iterator[tuple[str, object, str | None]]:
    for cls in PUBLIC_CLASSES:
        yield cls.__name__, cls, None
        for member_name, raw_member in vars(cls).items():
            if member_name.startswith("_"):
                continue
            member = raw_member
            if isinstance(member, (classmethod, staticmethod)):
                member = member.__func__
            elif isinstance(member, property):
                member = member.fget
            if not inspect.isfunction(member):
                continue
            yield (
                f"{cls.__name__}.{member_name}",
                member,
                f".{member_name}",
            )

    for label, function in PUBLIC_FUNCTIONS:
        yield label, function, f"{function.__name__}("


PUBLIC_DOCSTRING_TARGETS = tuple(iter_public_docstring_targets())


def test_public_docstring_target_matrix_is_complete() -> None:
    assert len(PUBLIC_CLASSES) == 27
    assert len(PUBLIC_FUNCTIONS) == 3
    assert len(PUBLIC_DOCSTRING_TARGETS) == EXPECTED_PUBLIC_TARGET_COUNT


def test_maintenance_helpers_are_not_public_api() -> None:
    assert not hasattr(TidyDataset, "build_manifest")
    assert not hasattr(ReactomeDb, "mapping_frame")
    assert not hasattr(WikiPathwaysDb, "lazy_frame")
    assert not hasattr(StringDb, "alias_schema")
    assert not hasattr(go_ontology, "build_tidy_frames")


# This is a structural floor. Whether an observed result explains why callers
# use the symbol still requires review against its call sites and producer tests.
@pytest.mark.parametrize(
    ("label", "target", "usage_token"),
    PUBLIC_DOCSTRING_TARGETS,
    ids=[target[0] for target in PUBLIC_DOCSTRING_TARGETS],
)
def test_public_api_docstring_has_direct_example(
    label: str,
    target: object,
    usage_token: str | None,
) -> None:
    docstring = inspect.getdoc(target)
    assert docstring is not None, f"{label} has no docstring"
    assert "Examples:" in docstring, f"{label} has no Examples section"

    examples = doctest.DocTestParser().get_examples(docstring)
    assert examples, f"{label} has no doctest prompts"
    assert any(example.want.strip() for example in examples), (
        f"{label} has no observable output or expected exception"
    )

    if usage_token is not None:
        sources = "\n".join(example.source for example in examples)
        assert usage_token in sources, (
            f"{label} does not demonstrate {usage_token!r} directly"
        )


@pytest.mark.parametrize(
    "target",
    EXECUTABLE_DOCSTRING_TARGETS,
    ids=[target.__name__ for target in EXECUTABLE_DOCSTRING_TARGETS],
)
def test_pure_public_examples_execute(target: type[object]) -> None:
    docstring = inspect.getdoc(target)
    assert docstring is not None
    module = sys.modules[target.__module__]
    test = doctest.DocTestParser().get_doctest(
        docstring,
        module.__dict__.copy(),
        target.__name__,
        inspect.getsourcefile(target) or target.__module__,
        0,
    )
    failures: list[str] = []
    result = doctest.DocTestRunner(
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    ).run(test, out=failures.append)
    assert result.failed == 0, "".join(failures)
