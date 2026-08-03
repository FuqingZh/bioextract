from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from .._shared import create_group_input_frames, create_input_id_frame
from .constant import (
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    OmniPathResourceName,
)
from .util import (
    extract_enzsub_frame,
    extract_interactions_frame,
)
from .util import (
    has_any_enzsub_modification as has_any_enzsub_modification_in_file,
)
from .util import (
    has_any_enzsub_relation as has_any_enzsub_relation_in_file,
)
from .util import (
    has_any_interaction_relation as has_any_interaction_relation_in_file,
)

__all__ = [
    "OmniPathDatabase",
]


@dataclass(frozen=True, slots=True)
class _OmniPathSnapshot:
    file_enzsub: Path | None = None
    file_interactions: Path | None = None


@dataclass(slots=True)
class OmniPathDatabase:
    """Path-first access to local OmniPath relation files.

    `OmniPathDatabase` is the public entrypoint for extracting OmniPath
    enzyme-substrate relations and interaction relations from local files.
    It exposes single and grouped selections through one selection type.

    Examples:
        Select enzyme-substrate relations for normalized protein IDs:

        >>> db = OmniPathDatabase.from_files(
        ...     enzsub="fixtures/omnipath/enzsub.tsv"
        ... )
        >>> (
        ...     db.select_ids(["P31749"])
        ...     .extract_enzsub()
        ...     .select("TargetId", "TargetSite")
        ...     .to_dicts()
        ... )
        [{'TargetId': 'BAD', 'TargetSite': 'S136'}, {'TargetId': 'FOXO3', 'TargetSite': 'T32'}]
    """

    snapshot: _OmniPathSnapshot

    @classmethod
    def from_files(
        cls,
        *,
        enzsub: os.PathLike[str] | str | None = None,
        interactions: os.PathLike[str] | str | None = None,
    ) -> OmniPathDatabase:
        """Create a dataset handle from local OmniPath files.

        Args:
            enzsub: Path to a local OmniPath `enzsub` text or gzip file.
            interactions: Path to a local OmniPath `interactions` text or
                gzip file.

        Returns:
            A dataset handle that can produce single or grouped selections.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If no resource files are provided.

        Examples:
            Open both relation resources from one fixture snapshot:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv",
            ...     interactions="fixtures/omnipath/interactions.tsv",
            ... )
            >>> sorted(db.available_resources)
            ['enzsub', 'interactions']
        """
        if enzsub is None and interactions is None:
            raise ValueError("At least one OmniPath resource file must be provided")

        file_enzsub = enzsub
        file_interactions = interactions
        if file_enzsub is not None:
            file_enzsub = Path(file_enzsub)
            if not file_enzsub.exists():
                raise FileNotFoundError(
                    f"OmniPath enzsub file not found: {file_enzsub}"
                )

        if file_interactions is not None:
            file_interactions = Path(file_interactions)
            if not file_interactions.exists():
                raise FileNotFoundError(
                    f"OmniPath interactions file not found: {file_interactions}"
                )

        return cls(
            snapshot=_OmniPathSnapshot(
                file_enzsub=file_enzsub,
                file_interactions=file_interactions,
            ),
        )

    @property
    def available_resources(self) -> frozenset[OmniPathResourceName]:
        """Return the relation resources backed by files in this snapshot.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> sorted(db.available_resources)
            ['enzsub']
        """
        resources: set[OmniPathResourceName] = set()
        if self.snapshot.file_enzsub is not None:
            resources.add("enzsub")
        if self.snapshot.file_interactions is not None:
            resources.add("interactions")
        return frozenset(resources)

    def has_any_enzsub_relation(self) -> bool:
        """Report whether the enzsub file contains a valid relation.

        Raises:
            ValueError: If this snapshot has no enzsub file.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> db.has_any_enzsub_relation()
            True
        """
        if self.snapshot.file_enzsub is None:
            raise ValueError("Cannot inspect OmniPath enzsub without enzsub file")
        return has_any_enzsub_relation_in_file(file_enzsub=self.snapshot.file_enzsub)

    def has_any_interaction_relation(self) -> bool:
        """Report whether the interactions file contains a valid relation.

        Raises:
            ValueError: If this snapshot has no interactions file.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     interactions="fixtures/omnipath/interactions.tsv"
            ... )
            >>> db.has_any_interaction_relation()
            True
        """
        if self.snapshot.file_interactions is None:
            raise ValueError(
                "Cannot inspect OmniPath interactions without interactions file"
            )
        return has_any_interaction_relation_in_file(
            file_interactions=self.snapshot.file_interactions
        )

    def has_any_enzsub_modification(self, modification: str) -> bool:
        """Report whether enzsub contains a relation for one modification.

        Args:
            modification: Exact OmniPath modification label to match after the
                file's standard string normalization.

        Raises:
            ValueError: If this snapshot has no enzsub file.

        Examples:
            Test present and absent modification labels:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> db.has_any_enzsub_modification("phosphorylation")
            True
            >>> db.has_any_enzsub_modification("ubiquitination")
            False
        """
        if self.snapshot.file_enzsub is None:
            raise ValueError("Cannot inspect OmniPath enzsub without enzsub file")
        return has_any_enzsub_modification_in_file(
            file_enzsub=self.snapshot.file_enzsub,
            modification=modification,
        )

    def select_ids(self, ids: Iterable[str]) -> OmniPathSelection:
        """Create a single-query selection from normalized protein IDs.

        The selection initially enables every resource backed by this snapshot.

        Args:
            ids: Protein identifiers. Whitespace, blank values, duplicates, and
                UniProt pipe-style IDs follow the shared normalization contract.

        Examples:
            Select the substrates linked to one kinase:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> db.select_ids(["P31749"]).extract_enzsub()["TargetId"].to_list()
            ['BAD', 'FOXO3']
        """
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        return OmniPathSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_group_membership=None,
            _df_groups=None,
            resources_selected=self.available_resources,
        )

    def select_groups(
        self,
        ids_by_group: Mapping[str, Iterable[str]],
    ) -> OmniPathSelection:
        """Create a grouped selection that preserves ``GroupId`` in outputs.

        The selection initially enables every resource backed by this snapshot.
        Relation endpoints are matched once against globally unique normalized
        IDs, then matched relations are expanded through group membership.

        Args:
            ids_by_group: Group labels mapped to protein identifiers. Labels
                must remain non-empty and unique after stripping whitespace.

        Examples:
            Preserve comparison labels in a grouped selection:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> (
            ...     db.select_groups({"up": ["P31749"], "down": ["MAPK1"]})
            ...     .extract_enzsub()
            ...     .select("GroupId", "SourceId")
            ...     .unique()
            ...     .sort("GroupId")
            ...     .to_dicts()
            ... )
            [{'GroupId': 'down', 'SourceId': 'MAPK1'}, {'GroupId': 'up', 'SourceId': 'P31749'}]
        """
        grp_in_frames = create_group_input_frames(
            ids_by_group,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        return OmniPathSelection(
            dataset=self,
            _df_groups=grp_in_frames.df_groups,
            _df_input_ids=grp_in_frames.df_input_ids,
            _df_group_membership=grp_in_frames.df_group_membership,
            resources_selected=self.available_resources,
        )


@dataclass(slots=True)
class OmniPathSelection:
    """Selection handle for both single and grouped OmniPath queries.

    Examples:
        Use a returned selection to materialize substrate relations:

        >>> db = OmniPathDatabase.from_files(
        ...     enzsub="fixtures/omnipath/enzsub.tsv"
        ... )
        >>> selection = db.select_ids(["P31749"])
        >>> selection.extract_enzsub()["TargetId"].to_list()
        ['BAD', 'FOXO3']
    """

    dataset: OmniPathDatabase
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_group_membership: pl.DataFrame | None = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    resources_selected: frozenset[OmniPathResourceName]
    _df_enzsub: pl.DataFrame | None = field(default=None, repr=False)
    _df_interactions: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> db.select_ids(["P31749"]).is_grouped
            False
        """
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the group ID column when this selection is grouped."""
        return ("GroupId",) if self.is_grouped else ()

    def with_resources(
        self,
        resources: Iterable[OmniPathResourceName],
    ) -> OmniPathSelection:
        """Create a selection constrained to enzsub and/or interactions.

        Existing cached frames are retained only for resources that remain
        selected. File availability is checked when extraction is requested.

        Args:
            resources: Non-empty iterable containing ``"enzsub"``,
                ``"interactions"``, or both.

        Returns:
            A new selection with the same normalized IDs and group mode.

        Raises:
            ValueError: If the iterable is empty or contains another name.

        Examples:
            Restrict unmapped-ID reporting to enzyme-substrate relations:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv",
            ...     interactions="fixtures/omnipath/interactions.tsv",
            ... )
            >>> selection = db.select_ids(["P31749", "ERBB2"])
            >>> (
            ...     selection.with_resources(["enzsub"])
            ...     .extract_unmatched_ids()
            ...     .to_dicts()
            ... )
            [{'InputId': 'ERBB2'}]
        """
        resources_selected = frozenset(resources)
        resources_invalid = resources_selected.difference({"enzsub", "interactions"})
        if resources_invalid:
            raise ValueError(
                f"Unsupported OmniPath resources: {sorted(resources_invalid)}"
            )
        if not resources_selected:
            raise ValueError("At least one OmniPath resource must be selected")

        return OmniPathSelection(
            dataset=self.dataset,
            _df_input_ids=self._df_input_ids,
            _df_group_membership=self._df_group_membership,
            _df_groups=self._df_groups,
            resources_selected=resources_selected,
            _df_enzsub=self._df_enzsub if "enzsub" in resources_selected else None,
            _df_interactions=(
                self._df_interactions if "interactions" in resources_selected else None
            ),
        )

    def with_enzsub(self) -> OmniPathSelection:
        """Create a new selection constrained to OmniPath enzsub relations.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv",
            ...     interactions="fixtures/omnipath/interactions.tsv",
            ... )
            >>> (
            ...     db.select_ids(["P31749"])
            ...     .with_enzsub()
            ...     .extract_enzsub()["TargetId"]
            ...     .to_list()
            ... )
            ['BAD', 'FOXO3']
        """
        return self.with_resources(["enzsub"])

    def with_interactions(self) -> OmniPathSelection:
        """Create a new selection constrained to OmniPath interactions.

        Examples:
            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv",
            ...     interactions="fixtures/omnipath/interactions.tsv",
            ... )
            >>> (
            ...     db.select_ids(["ERBB2"])
            ...     .with_interactions()
            ...     .extract_interactions()
            ...     .select("SourceId", "TargetId")
            ...     .to_dicts()
            ... )
            [{'SourceId': 'EGFR', 'TargetId': 'ERBB2'}]
        """
        return self.with_resources(["interactions"])

    def extract_enzsub(self) -> pl.DataFrame:
        """Extract matched OmniPath enzyme-substrate relations.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `SourceId`, `TargetId`, `TargetSite`, `Modification`
            - grouped selection: `GroupId`, `SourceId`, `TargetId`, `TargetSite`,
              `Modification`

        Raises:
            ValueError: If the selection does not enable `enzsub`, if the
                `enzsub` file is missing, or if the file is missing required
                columns.

        Examples:
            Extract a substrate site and modification:

            >>> db = OmniPathDatabase.from_files(
            ...     enzsub="fixtures/omnipath/enzsub.tsv"
            ... )
            >>> db.select_ids(["P31749"]).extract_enzsub().head(1).to_dicts()
            [{'SourceId': 'P31749', 'TargetId': 'BAD', 'TargetSite': 'S136', 'Modification': 'phosphorylation'}]
        """
        if "enzsub" not in self.resources_selected:
            raise ValueError(
                "OmniPath resource 'enzsub' is not enabled for this selection"
            )
        if self.dataset.snapshot.file_enzsub is None:
            raise ValueError("Cannot extract OmniPath enzsub without enzsub file")
        if self._df_enzsub is None:
            self._df_enzsub = extract_enzsub_frame(
                file_enzsub=self.dataset.snapshot.file_enzsub,
                df_input_ids=self._df_input_ids,
                df_group_membership=self._df_group_membership,
            )
        return self._df_enzsub

    def extract_interactions(self) -> pl.DataFrame:
        """Extract matched OmniPath interaction relations.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `SourceId`, `TargetId`, `IsDirected`,
              `IsStimulation`, `IsInhibition`
            - grouped selection: `GroupId`, `SourceId`, `TargetId`,
              `IsDirected`, `IsStimulation`, `IsInhibition`

        Raises:
            ValueError: If the selection does not enable `interactions`, if the
                interactions file is missing, or if the file is missing
                required columns.

        Examples:
            Extract one interaction and its direction flags:

            >>> db = OmniPathDatabase.from_files(
            ...     interactions="fixtures/omnipath/interactions.tsv"
            ... )
            >>> db.select_ids(["ERBB2"]).extract_interactions().to_dicts()
            [{'SourceId': 'EGFR', 'TargetId': 'ERBB2', 'IsDirected': False, 'IsStimulation': True, 'IsInhibition': False}]
        """
        if "interactions" not in self.resources_selected:
            raise ValueError(
                "OmniPath resource 'interactions' is not enabled for this selection"
            )
        if self.dataset.snapshot.file_interactions is None:
            raise ValueError(
                "Cannot extract OmniPath interactions without interactions file"
            )
        if self._df_interactions is None:
            self._df_interactions = extract_interactions_frame(
                file_interactions=self.dataset.snapshot.file_interactions,
                df_input_ids=self._df_input_ids,
                df_group_membership=self._df_group_membership,
            )
        return self._df_interactions

    def extract_unmatched_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs not found in the selected resources.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `InputId`
            - grouped selection: `GroupId`, `InputId`

        Examples:
            Report an identifier absent from the selected resource:

            >>> db = OmniPathDatabase.from_files(
            ...     interactions="fixtures/omnipath/interactions.tsv"
            ... )
            >>> selection = db.select_ids(["MISSING"]).with_interactions()
            >>> selection.extract_unmatched_ids().to_dicts()
            [{'InputId': 'MISSING'}]
        """
        if self._df_unmapped is None:
            col_group_id = list(self._col_group_id)
            cols_index = col_group_id + ["InputId"]
            df_input_rows = self._df_input_ids
            if self.is_grouped:
                if self._df_group_membership is None:
                    raise RuntimeError(
                        "Grouped OmniPath selection lacks group membership"
                    )
                df_input_rows = self._df_group_membership

            df_matched_parts: list[pl.DataFrame] = []
            if "enzsub" in self.resources_selected:
                df_enzsub = self.extract_enzsub()
                df_matched_parts.extend(
                    [
                        df_enzsub.select(
                            col_group_id + [pl.col("SourceId").alias("InputId")]
                        ),
                        df_enzsub.select(
                            col_group_id + [pl.col("TargetId").alias("InputId")]
                        ),
                    ]
                )

            if "interactions" in self.resources_selected:
                df_interactions = self.extract_interactions()
                df_matched_parts.extend(
                    [
                        df_interactions.select(
                            col_group_id + [pl.col("SourceId").alias("InputId")]
                        ),
                        df_interactions.select(
                            col_group_id + [pl.col("TargetId").alias("InputId")]
                        ),
                    ]
                )

            df_matched_input_ids = (
                pl.concat(df_matched_parts, how="vertical_relaxed")
                .join(df_input_rows, on=cols_index, how="inner")
                .unique(subset=cols_index)
                .sort(cols_index)
            )
            self._df_unmapped = (
                df_input_rows.join(
                    df_matched_input_ids,
                    on=cols_index,
                    how="anti",
                )
                .select(cols_index)
                .sort(cols_index)
            )
        return self._df_unmapped
