from __future__ import annotations

import doctest
import inspect
from collections.abc import Iterator
from types import FunctionType
from typing import Protocol, cast

import pytest

from bioextract._tidy import TidyDataset
from bioextract.chebi import ChEBIDatabase
from bioextract.eggnog import EggNOGDatabase, EggnogSelection
from bioextract.go import GODatabase, GoSubsetId
from bioextract.interpro import (
    InterProDatabase,
    InterProSelection,
)
from bioextract.kegg import KEGGDatabase
from bioextract.kegg.brite import (
    build_tidy_frames as build_kegg_tidy_frames,
)
from bioextract.kegg.kegg import KeggSelection
from bioextract.omnipath import OmniPathDatabase
from bioextract.omnipath.omnipath import OmniPathSelection
from bioextract.reactome import ReactomeDatabase
from bioextract.reactome.reactome import ReactomeSelection
from bioextract.rhea import RheaDatabase, RheaReactionSelection, RheaWriteResult
from bioextract.stringdb import STRINGDatabase
from bioextract.stringdb.stringdb import StringSelection
from bioextract.uniprot import UniProtDatabase
from bioextract.wikipathways import WikiPathwaysDatabase
from bioextract.wikipathways.wikipathways import WikiPathwaysSelection

PUBLIC_CLASSES = (
    ChEBIDatabase,
    GODatabase,
    GoSubsetId,
    KEGGDatabase,
    KeggSelection,
    ReactomeDatabase,
    ReactomeSelection,
    WikiPathwaysDatabase,
    WikiPathwaysSelection,
    EggNOGDatabase,
    EggnogSelection,
    InterProDatabase,
    InterProSelection,
    UniProtDatabase,
    STRINGDatabase,
    StringSelection,
    OmniPathDatabase,
    OmniPathSelection,
    RheaDatabase,
    RheaReactionSelection,
    RheaWriteResult,
    TidyDataset,
)

PUBLIC_FUNCTIONS = (("kegg.brite.build_tidy_frames", build_kegg_tidy_frames),)

EXPECTED_PUBLIC_TARGET_COUNT = 148


class _FunctionDescriptor(Protocol):
    __func__: FunctionType


class _FunctionProperty(Protocol):
    fget: FunctionType | None


def iter_public_docstring_targets() -> Iterator[tuple[str, object, str | None]]:
    for cls in PUBLIC_CLASSES:
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

    for label, function in PUBLIC_FUNCTIONS:
        yield label, function, f"{function.__name__}("


PUBLIC_DOCSTRING_TARGETS = tuple(iter_public_docstring_targets())


def test_public_docstring_target_matrix_is_complete() -> None:
    assert len(PUBLIC_CLASSES) == 22
    assert len(PUBLIC_FUNCTIONS) == 1
    assert len(PUBLIC_DOCSTRING_TARGETS) == EXPECTED_PUBLIC_TARGET_COUNT


def test_maintenance_helpers_are_not_public_api() -> None:
    assert not hasattr(TidyDataset, "build_manifest")
    assert not hasattr(ReactomeDatabase, "mapping_frame")
    assert not hasattr(WikiPathwaysDatabase, "lazy_frame")
    assert not hasattr(STRINGDatabase, "alias_schema")


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
