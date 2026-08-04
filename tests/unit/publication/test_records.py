from dataclasses import FrozenInstanceError, fields

import pytest

from bioextract.publication import (
    PublicationColumnMapping,
    PublicationInspection,
    PublicationMetadata,
    PublicationSourceFile,
    PublicationTable,
    PublicationValidationIssue,
)


def test_publication_records_are_immutable_and_slotted() -> None:
    records = (
        PublicationMetadata("key", "value"),
        PublicationTable("term", "canonical", 1),
        PublicationSourceFile("source", "input.tsv", 10, "text/tab-separated", None),
        PublicationColumnMapping("term", "Term ID", "term_id", "normalized"),
        PublicationValidationIssue(
            1,
            "warning",
            "fixture",
            "source",
            "term",
            None,
            None,
            None,
            None,
            None,
            "fixture issue",
        ),
    )
    for record in records:
        assert not hasattr(record, "__dict__")
        with pytest.raises(FrozenInstanceError):
            record.__setattr__(fields(record)[0].name, "changed")


def test_publication_inspection_has_exact_stable_fields() -> None:
    assert tuple(field.name for field in fields(PublicationInspection)) == (
        "path",
        "resource_name",
        "resource_schema_version",
        "source_schema_profile",
        "source_schema_version",
        "release_version",
        "release_version_source",
        "metadata_schema_version",
        "package_version",
        "generated_at",
        "scope",
        "validation_status",
        "validation_issue_count",
        "table_counts_verified",
        "metadata",
        "tables",
        "source_files",
        "column_mappings",
        "validation_issues",
    )
    assert "__dict__" not in PublicationInspection.__slots__
