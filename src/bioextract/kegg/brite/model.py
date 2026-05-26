from dataclasses import dataclass, field

import polars as pl


@dataclass(slots=True)
class PathwayLevelRecord:
    id: str
    name: str | None
    kegg_id: str | None = None


@dataclass(slots=True)
class PathwayLeafRecord:
    entry: PathwayLevelRecord
    ko: PathwayLevelRecord | None


@dataclass(slots=True)
class BriteRecord:
    pathway_level1_id: str
    pathway_level1_name: str | None
    pathway_level2_id: str
    pathway_level2_name: str | None
    pathway_level3_id: str
    pathway_level3_kegg_id: str | None
    pathway_level3_name: str | None
    entry_id: str | None
    entry_name: str | None
    ko_id: str | None
    ko_name: str | None


@dataclass(slots=True)
class BriteColumnBuffer:
    pathway_level1_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level1_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level2_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level2_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level3_id: list[str] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level3_kegg_id: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    pathway_level3_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    entry_id: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    entry_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    ko_id: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )
    ko_name: list[str | None] = field(
        default_factory=lambda: [], metadata={"dtype": pl.String}
    )

    def append_record(self, record: BriteRecord) -> None:
        self.pathway_level1_id.append(record.pathway_level1_id)
        self.pathway_level1_name.append(record.pathway_level1_name)
        self.pathway_level2_id.append(record.pathway_level2_id)
        self.pathway_level2_name.append(record.pathway_level2_name)
        self.pathway_level3_id.append(record.pathway_level3_id)
        self.pathway_level3_kegg_id.append(record.pathway_level3_kegg_id)
        self.pathway_level3_name.append(record.pathway_level3_name)
        self.entry_id.append(record.entry_id)
        self.entry_name.append(record.entry_name)
        self.ko_id.append(record.ko_id)
        self.ko_name.append(record.ko_name)
