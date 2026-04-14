import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from bioextract._shared import validate_count_limit, validate_file_size

from .._shared import create_group_input_frames, create_input_id_frame
from .constant import (
    SCHEMA_GROUP_INPUT_IDS,
    SCHEMA_GROUPS,
    SCHEMA_UNMAPPED,
    OmniPathResourceName,
)
from .spec import OmniPathResourceLimits
from .util import (
    extract_enzsub_frame,
    extract_interactions_frame,
    has_any_enzsub_modification as has_any_enzsub_modification_in_file,
    has_any_enzsub_relation as has_any_enzsub_relation_in_file,
    has_any_interaction_relation as has_any_interaction_relation_in_file,
)

__all__ = [
    "OmniPathDb",
]


@dataclass(frozen=True, slots=True)
class _OmniPathSnapshot:
    file_enzsub: Path | None = None
    file_interactions: Path | None = None


@dataclass(slots=True)
class OmniPathDb:
    """Path-first access to local OmniPath relation files.

    `OmniPathDb` is the public entrypoint for extracting OmniPath
    enzyme-substrate relations and interaction relations from local files.
    It keeps dataset-level resource limits and exposes single and grouped
    selections through one selection type.
    """

    snapshot: _OmniPathSnapshot
    limits: OmniPathResourceLimits = field(default_factory=OmniPathResourceLimits)

    DEFAULT_RESOURCE_LIMITS = OmniPathResourceLimits()

    @classmethod
    def from_files(
        cls,
        *,
        file_enzsub: os.PathLike[str] | str | None = None,
        file_interactions: os.PathLike[str] | str | None = None,
        limits: OmniPathResourceLimits | None = None,
    ) -> "OmniPathDb":
        """Create a dataset handle from local OmniPath files.

        Args:
            file_enzsub: Path to a local OmniPath `enzsub` text or gzip file.
            file_interactions: Path to a local OmniPath `interactions` text or
                gzip file.
            limits: Dataset-level resource limits. When omitted, default
                fail-fast limits are used. see :class:`OmniPathResourceLimits` and
                `OmniPathDb.DEFAULT_RESOURCE_LIMITS` for details.

        Returns:
            A dataset handle that can produce single or grouped selections.

        Raises:
            FileNotFoundError: If any provided file does not exist.
            ValueError: If no resource files are provided or a configured
                file-size limit is exceeded.
        """
        if file_enzsub is None and file_interactions is None:
            raise ValueError("At least one OmniPath resource file must be provided")

        limits_resolved = OmniPathResourceLimits() if limits is None else limits

        if file_enzsub is not None:
            file_enzsub = Path(file_enzsub)
            if not file_enzsub.exists():
                raise FileNotFoundError(
                    f"OmniPath enzsub file not found: {file_enzsub}"
                )
            validate_file_size(
                file_path=file_enzsub,
                size_max=limits_resolved.file_enzsub_bytes_max,
                label="OmniPath enzsub file",
            )

        if file_interactions is not None:
            file_interactions = Path(file_interactions)
            if not file_interactions.exists():
                raise FileNotFoundError(
                    f"OmniPath interactions file not found: {file_interactions}"
                )
            validate_file_size(
                file_path=file_interactions,
                size_max=limits_resolved.file_interactions_bytes_max,
                label="OmniPath interactions file",
            )

        return cls(
            snapshot=_OmniPathSnapshot(
                file_enzsub=file_enzsub,
                file_interactions=file_interactions,
            ),
            limits=limits_resolved,
        )

    @property
    def available_resources(self) -> frozenset[OmniPathResourceName]:
        resources: set[OmniPathResourceName] = set()
        if self.snapshot.file_enzsub is not None:
            resources.add("enzsub")
        if self.snapshot.file_interactions is not None:
            resources.add("interactions")
        return frozenset(resources)

    def has_any_enzsub_relation(self) -> bool:
        """Report whether the dataset contains any valid enzsub relation."""
        if self.snapshot.file_enzsub is None:
            raise ValueError("Cannot inspect OmniPath enzsub without enzsub file")
        return has_any_enzsub_relation_in_file(file_enzsub=self.snapshot.file_enzsub)

    def has_any_interaction_relation(self) -> bool:
        """Report whether the dataset contains any valid interaction relation."""
        if self.snapshot.file_interactions is None:
            raise ValueError(
                "Cannot inspect OmniPath interactions without interactions file"
            )
        return has_any_interaction_relation_in_file(
            file_interactions=self.snapshot.file_interactions
        )

    def has_any_enzsub_modification(self, modification: str) -> bool:
        """Report whether the dataset contains any valid enzsub relation for a modification."""
        if self.snapshot.file_enzsub is None:
            raise ValueError("Cannot inspect OmniPath enzsub without enzsub file")
        return has_any_enzsub_modification_in_file(
            file_enzsub=self.snapshot.file_enzsub,
            modification=modification,
        )

    def select_ids(self, ids: Iterable[str]) -> OmniPathSelection:
        """Create a single-query selection from input IDs."""
        df_input_ids = create_input_id_frame(ids, schema_unmapped=SCHEMA_UNMAPPED)
        validate_count_limit(
            count=df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return OmniPathSelection(
            dataset=self,
            _df_input_ids=df_input_ids,
            _df_groups=None,
            resources_selected=self.available_resources,
        )

    def select_groups(
        self,
        group_to_ids: Mapping[str, Iterable[str]],
    ) -> OmniPathSelection:
        """Create a grouped selection from multiple input-ID sets."""
        grp_in_frames = create_group_input_frames(
            group_to_ids,
            schema_groups=SCHEMA_GROUPS,
            schema_group_input_ids=SCHEMA_GROUP_INPUT_IDS,
        )
        validate_count_limit(
            count=grp_in_frames.df_groups.height,
            limit_max=self.limits.num_groups_max,
            label="Group count",
        )
        validate_count_limit(
            count=grp_in_frames.df_input_ids.height,
            limit_max=self.limits.num_input_ids_max,
            label="Normalized input ID count",
        )
        return OmniPathSelection(
            dataset=self,
            _df_groups=grp_in_frames.df_groups,
            _df_input_ids=grp_in_frames.df_input_ids,
            resources_selected=self.available_resources,
        )


@dataclass(slots=True)
class OmniPathSelection:
    """Selection handle for both single and grouped OmniPath queries."""

    dataset: OmniPathDb
    _df_input_ids: pl.DataFrame = field(repr=False)
    _df_groups: pl.DataFrame | None = field(repr=False)
    resources_selected: frozenset[OmniPathResourceName]
    _df_enzsub: pl.DataFrame | None = field(default=None, repr=False)
    _df_interactions: pl.DataFrame | None = field(default=None, repr=False)
    _df_unmapped: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def is_grouped(self) -> bool:
        """Report whether this selection carries `GroupId` through outputs."""
        return self._df_groups is not None

    @property
    def _col_group_id(self) -> tuple[str, ...]:
        """Return the group ID column when this selection is grouped."""
        return ("GroupId",) if self.is_grouped else ()

    def with_resources(
        self,
        resources: Iterable[OmniPathResourceName],
    ) -> OmniPathSelection:
        """Create a new selection constrained to the given OmniPath resources."""
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
            _df_groups=self._df_groups,
            resources_selected=resources_selected,
            _df_enzsub=self._df_enzsub if "enzsub" in resources_selected else None,
            _df_interactions=(
                self._df_interactions if "interactions" in resources_selected else None
            ),
        )

    def with_enzsub(self) -> OmniPathSelection:
        """Create a new selection constrained to OmniPath enzsub relations."""
        return self.with_resources(["enzsub"])

    def with_interactions(self) -> OmniPathSelection:
        """Create a new selection constrained to OmniPath interactions."""
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
                cols_group_id=self._col_group_id,
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
                cols_group_id=self._col_group_id,
            )
        return self._df_interactions

    def extract_unmapped_input_ids(self) -> pl.DataFrame:
        """Extract normalized input IDs not found in the selected resources.

        Returns:
            A materialized table with one of these schemas:

            - single selection: `InputId`
            - grouped selection: `GroupId`, `InputId`
        """
        if self._df_unmapped is None:
            col_group_id = list(self._col_group_id)
            cols_index = col_group_id + ["InputId"]

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
                .join(self._df_input_ids, on=cols_index, how="inner")
                .unique(subset=cols_index)
                .sort(cols_index)
            )
            self._df_unmapped = (
                self._df_input_ids.join(
                    df_matched_input_ids,
                    on=cols_index,
                    how="anti",
                )
                .select(cols_index)
                .sort(cols_index)
            )
        return self._df_unmapped
