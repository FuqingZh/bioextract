from __future__ import annotations

import re
from pathlib import Path

import pytest

import bioextract
from bioextract.chebi._query import ChEBICompoundSelection
from bioextract.eggnog.eggnog import EggnogSelection
from bioextract.go._query import GOAncestorSelection
from bioextract.kegg.metabolic.core import KEGGMetabolicSelection
from bioextract.omnipath.omnipath import OmniPathSelection
from bioextract.rhea._query import RheaReactionSelection
from bioextract.stringdb.stringdb import StringSelection
from bioextract.uniprot._query import UniProtSelection
from bioextract.uniprot.uniprot import UniProtDatabase

REPOSITORY_ROOT = Path(__file__).parents[3]
README_PATHS = (
    REPOSITORY_ROOT / "README.md",
    REPOSITORY_ROOT / "README.zh-CN.md",
)
PYTHON_FENCE = re.compile(r"```python\n(?P<source>.*?)\n```", re.DOTALL)
MARKDOWN_LINK = re.compile(r"\[[^]]+\]\((?P<target>[^)]+)\)")

PUBLIC_OWNER_NAMES = (
    "ChEBIDatabase",
    "EggNOGDatabase",
    "GODatabase",
    "InterProDatabase",
    "KEGGDatabase",
    "OmniPathDatabase",
    "ReactomeDatabase",
    "RheaDatabase",
    "STRINGDatabase",
    "UniProtDatabase",
    "WikiPathwaysDatabase",
    "inspect_publication",
)

EXPECTED_SNIPPETS = (
    "selection.compounds()",
    "selection.names()",
    "selection.relations()",
    "selection.matches()",
    "selection.reactions()",
    "selection.participants()",
    "selection.cross_references()",
    "selection.ancestors()",
    "selection.pathway_memberships()",
    'db.select_ids(["9606.ENSP00000369497"]).mappings()',
    'mapping.scan_mapping(taxon_ids=["9606"])',
    ").proteins()",
    "selection.mappings()",
    "selection.edges()",
    "selection.enzsub()",
    "selection.unmatched_ids()",
)

TERMINAL_OWNERS = (
    (ChEBICompoundSelection, ("compounds", "names", "relations", "unmatched_ids")),
    (
        RheaReactionSelection,
        (
            "matches",
            "reactions",
            "participants",
            "cross_references",
            "unmatched_ids",
        ),
    ),
    (GOAncestorSelection, ("ancestors", "unmatched_ids")),
    (
        KEGGMetabolicSelection,
        ("reactions", "pathway_memberships", "unmatched_ids"),
    ),
    (EggnogSelection, ("mappings",)),
    (UniProtDatabase, ("scan_mapping",)),
    (UniProtSelection, ("proteins",)),
    (StringSelection, ("mappings", "unmatched_ids", "edges")),
    (OmniPathSelection, ("enzsub", "unmatched_ids")),
)


@pytest.mark.parametrize("readme_path", README_PATHS, ids=lambda path: path.name)
def test_readme_python_examples_use_current_api(readme_path: Path) -> None:
    text = readme_path.read_text(encoding="utf-8")
    assert ".extract_" not in text
    assert ".read_mapping(" not in text
    for snippet in EXPECTED_SNIPPETS:
        assert snippet in text

    sources = [match.group("source") for match in PYTHON_FENCE.finditer(text)]
    assert sources
    for index, source in enumerate(sources, start=1):
        compile(source, f"{readme_path.name}:python-block-{index}", "exec")


@pytest.mark.parametrize("readme_path", README_PATHS, ids=lambda path: path.name)
def test_readme_local_links_resolve(readme_path: Path) -> None:
    text = readme_path.read_text(encoding="utf-8")
    for match in MARKDOWN_LINK.finditer(text):
        target = match.group("target")
        if target.startswith(("#", "http://", "https://", "mailto:")):
            continue
        path_text = target.split("#", maxsplit=1)[0]
        assert (readme_path.parent / path_text).exists(), target


@pytest.mark.parametrize(
    ("owner", "terminal_names"),
    TERMINAL_OWNERS,
    ids=lambda value: value.__name__ if isinstance(value, type) else None,
)
def test_readme_noun_terminals_are_importable(
    owner: type[object],
    terminal_names: tuple[str, ...],
) -> None:
    for terminal_name in terminal_names:
        assert callable(getattr(owner, terminal_name, None))


@pytest.mark.parametrize("owner_name", PUBLIC_OWNER_NAMES)
def test_readme_public_owners_resolve_from_package(owner_name: str) -> None:
    assert callable(getattr(bioextract, owner_name, None))
