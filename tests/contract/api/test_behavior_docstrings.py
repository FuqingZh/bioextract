from __future__ import annotations

import doctest
import inspect
from collections.abc import Iterator
from types import FunctionType
from typing import Protocol, cast

import pytest

from bioextract import (
    ChEBIDatabase,
    EggNOGDatabase,
    GODatabase,
    InterProDatabase,
    KEGGDatabase,
    OmniPathDatabase,
    ReactomeDatabase,
    RheaDatabase,
    STRINGDatabase,
    UniProtDatabase,
    WikiPathwaysDatabase,
)
from bioextract._tidy import TidyDataset
from bioextract.chebi._query import ChEBICompoundSelection
from bioextract.eggnog.eggnog import EggnogSelection
from bioextract.go.go import GoSubsetId
from bioextract.interpro.interpro import InterProSelection
from bioextract.kegg.brite.tidy import (
    build_tidy_frames as build_kegg_tidy_frames,
)
from bioextract.kegg.kegg import KeggSelection
from bioextract.kegg.metabolic.core import KEGGMetabolicSelection
from bioextract.omnipath.omnipath import OmniPathSelection
from bioextract.reactome.reactome import ReactomeSelection
from bioextract.rhea._query import RheaReactionSelection
from bioextract.rhea.rhea import RheaWriteResult
from bioextract.stringdb.stringdb import StringSelection
from bioextract.uniprot._query import UniProtSelection
from bioextract.wikipathways.wikipathways import WikiPathwaysSelection

BEHAVIOR_CLASSES = (
    ChEBIDatabase,
    ChEBICompoundSelection,
    GODatabase,
    GoSubsetId,
    KEGGDatabase,
    KeggSelection,
    KEGGMetabolicSelection,
    ReactomeDatabase,
    ReactomeSelection,
    WikiPathwaysDatabase,
    WikiPathwaysSelection,
    EggNOGDatabase,
    EggnogSelection,
    InterProDatabase,
    InterProSelection,
    UniProtDatabase,
    UniProtSelection,
    STRINGDatabase,
    StringSelection,
    OmniPathDatabase,
    OmniPathSelection,
    RheaDatabase,
    RheaReactionSelection,
    RheaWriteResult,
    TidyDataset,
)

BEHAVIOR_FUNCTIONS = (("kegg.brite.build_tidy_frames", build_kegg_tidy_frames),)

EXPECTED_BEHAVIOR_TARGET_COUNT = 186


class _FunctionDescriptor(Protocol):
    __func__: FunctionType


class _FunctionProperty(Protocol):
    fget: FunctionType | None


def iter_behavior_docstring_targets() -> Iterator[tuple[str, object, str | None]]:
    for cls in BEHAVIOR_CLASSES:
        yield cls.__name__, cls, None
        for member_name, raw_member in vars(cls).items():
            if member_name.startswith("_"):
                continue
            member = raw_member
            if isinstance(member, (classmethod, staticmethod)):
                member = cast(_FunctionDescriptor, member).__func__
            elif isinstance(member, property):
                member = cast(_FunctionProperty, member).fget
            if not inspect.isfunction(member):
                continue
            yield (
                f"{cls.__name__}.{member_name}",
                member,
                f".{member_name}",
            )

    for label, function in BEHAVIOR_FUNCTIONS:
        yield label, function, f"{function.__name__}("


BEHAVIOR_DOCSTRING_TARGETS = tuple(iter_behavior_docstring_targets())


def test_behavior_docstring_target_matrix_is_complete() -> None:
    assert len(BEHAVIOR_CLASSES) == 25
    assert len(BEHAVIOR_FUNCTIONS) == 1
    assert len(BEHAVIOR_DOCSTRING_TARGETS) == EXPECTED_BEHAVIOR_TARGET_COUNT


def test_maintenance_helpers_are_not_public_api() -> None:
    assert not hasattr(TidyDataset, "build_manifest")
    assert not hasattr(ReactomeDatabase, "mapping_frame")
    assert not hasattr(WikiPathwaysDatabase, "lazy_frame")
    assert not hasattr(STRINGDatabase, "alias_schema")


# This is a structural floor. Whether an observed result explains why callers
# use the symbol still requires review against its call sites and producer tests.
@pytest.mark.parametrize(
    ("label", "target", "usage_token"),
    BEHAVIOR_DOCSTRING_TARGETS,
    ids=[target[0] for target in BEHAVIOR_DOCSTRING_TARGETS],
)
def test_behavior_docstring_has_direct_example(
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
