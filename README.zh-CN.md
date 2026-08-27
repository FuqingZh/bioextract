<p align="center">
  <a href="https://github.com/FuqingZh/bioextract/blob/main/README.md">English</a> |
  <strong>简体中文</strong>
</p>

# bioextract

面向官方生物数据库快照的稳定、来源可追溯的领域访问层。

`bioextract` 封装不同资源特有的文件布局、标识符规则、层级关系、方向语义和重复
连接操作。调用方提供本地快照文件；本库既不负责下载资源，也不感知最终由哪个
应用程序消费这些数据。

> 本文档是英文 README 的简体中文版本。若翻译与英文原文存在差异，以
> [英文版](README.md)为准。

## 架构

领域契约优先，存储只是执行策略：

- 当上游原生表示适合所支持的查询时，直接访问官方/原生数据（eggNOG、STRING
  和 OmniPath）；
- 每个需要物化的逻辑产品发布为一个由 bioextract 管理的 DuckDB，无论其中包含
  一个还是多个关系；
- 常规过滤、排序、分组和 SQL 留给 Polars 或 DuckDB；
- 仅当便捷方法承载资源自身定义的 ID 解析、关系遍历、分组或未匹配 ID 统计时，
  才将其加入领域 API。

物化写入器通过 `write_duckdb(path)` 接收明确的目标路径。读取器使用
`XDatabase.from_duckdb(path)`，每次调用 `connect()` 都会返回一个全新、由调用方
管理的只读 DuckDB 连接。写入器先在暂存文件中完成验证，只有成功后才原子发布。
发布来源信息恰好存放在 DuckDB `_bioextract` schema 的五个关系中；生物学关系
存放在 `main`。Metadata v2 是目前唯一受支持的发布元数据契约。Parquet 只作为
上游格式、内部传输格式或通用交换格式使用，绝不是 bioextract 的规范发布格式。

在添加资源、公共查询方法或存储策略之前，请先阅读
[领域访问架构](docs/architecture/20260729-v1.0-domain-access-architecture.md)。

## 检查发布产物

检查一个明确指定的本地发布产物，无需选择资源专用读取器，也无需扫描生物学
数据行：

```python
import bioextract

publication = bioextract.inspect_publication("out/go.duckdb")

print(publication.resource_name)
print(publication.resource_schema_version)
print(publication.release_version)
print(publication.validation_status)
print(publication.table_counts_verified)  # False by default
for table in publication.tables:
    print(table.table_name, table.table_role, table.row_count)
```

`inspect_publication()` 是稳定的顶层函数。不可变结果及其辅助记录类型仍可从
`bioextract.publication` 获取，但不作为包的顶层导出。

## 安装

```console
pip install bioextract
```

## ChEBI 与 ChemOnt

FULL OBO 提供规范的化合物、标识符、名称、交叉引用、属性和关系 schema。SDF
只补充 molfile 记录；可选的 ChemOnt 在同一容器中仍保持为独立的
`chemont_*` 图：

```python
from bioextract import ChEBIDatabase

result = ChEBIDatabase.from_obo(
    "chebi/database/2026-07-07/raw/chebi.obo",
    chemont_obo="ChemOnt_2_1.obo.zip",
).write_duckdb("out/chebi.duckdb")

print(result.tables)
```

打开发布产物后，可进行稳定的领域提取，或执行不受限制的原生只读 SQL：

```python
database = ChEBIDatabase.from_duckdb("out/chebi.duckdb")
selection = database.select_compounds(
    ["CHEBI:15377", "CHEBI:10743"],
    namespace="chebi",
)

lf_compounds = selection.compounds()
lf_names = selection.names()
lf_relations = selection.relations()
lf_unmatched = selection.unmatched_ids()

with database.connect() as connection:
    prefixes = connection.execute(
        "SELECT DISTINCT source_prefix FROM compound_cross_reference"
    ).fetchall()
```

外部交叉引用直接使用官方前缀作为 `namespace`，例如 `kegg.compound` 或 `hmdb`。
公开共享 ID 采用完整的 `CHEBI:<number>` CURIE。只有构建不完整来源时才使用显式
TSV 文件；适用时，本库会在内部识别纯文本、gzip、zip 和 tar 输入。

## Rhea

从完整解压的发布版本或归档文件构建一个可直接查询的数据库：

```python
from bioextract import RheaDatabase

result = RheaDatabase.from_files("rhea-release.zip").write_duckdb(
    "out/rhea.duckdb"
)
print(result.tables)
```

显式文件输入允许能力不完整或混合的组合，同时仍写入同一种 DuckDB 容器：

```python
from bioextract import RheaDatabase

RheaDatabase.from_files(
    rdf="rhea.rdf.gz",
    directions="rhea-directions.tsv",
    relationships="rhea-relationships.tsv",
    xrefs="rhea2xrefs.tsv",
    uniprot_sprot="rhea2uniprot_sprot.tsv",
    uniprot_trembl="rhea2uniprot_trembl.tsv.gz",
).write_duckdb("out/rhea.duckdb")
```

打开已发布数据库，并通过任一受支持的官方命名空间选择反应：

```python
database = RheaDatabase.from_duckdb("out/rhea.duckdb")
selection = database.select_reactions(
    ["CHEBI:15377", "CHEBI:16474"],
    namespace="chebi",
)

lf_matches = selection.matches()
lf_reactions = selection.reactions()
lf_participants = selection.participants()
lf_cross_references = selection.cross_references()
lf_unmatched = selection.unmatched_ids()
```

`select_reactions()` 和 `select_groups()` 是延迟执行的领域查询计划；其名词终端
返回可重放的 Polars `LazyFrame`。在需要立即求值 `DataFrame` 的应用边界调用
`.collect()`。参与物输出会保留准确的
Rhea ID、主反应 ID、方向、反应侧和化合物字段。ChEBI 字段使用完整的
`CHEBI:<number>` CURIE，无需拼接前缀或类型转换即可与 ChEBI 发布产物等值连接。
`directional_role` 只为 `LR` 和 `RL` 填充值；方向未定义或双向的反应保留空值，
而不会虚构底物/产物方向。

方向、层级、数据表和来源契约详见
[Rhea 架构](docs/architecture/rhea-db.md)。

## GO

GO 是多关系本体，整体发布为一个 DuckDB：

```python
from bioextract import GODatabase

go = GODatabase.from_obo("go-basic.obo")
df_terms = go.select_terms(subset_id="goslim_generic")
df_cellular_components = go.select_terms(namespace="cellular_component")
selection = go.select_ancestors(
    ["GO:0008150", "GO:1234567"],
    target_subset_id="goslim_generic",
    include_self=True,
)
lf_ancestors = selection.ancestors()
lf_unmatched = selection.unmatched_ids()
result = go.write_duckdb("out/go.duckdb")
```

数据表包括 `term`、`term_relation`、`term_synonym`、`term_xref`、
`term_alternate_id`、`term_ancestor` 和 `term_depth`。

GO 祖先选择会解析规范或替代 GO ID，并可将其 `is_a`/`part_of` 祖先投影到一个
OBO 子集。蛋白质成员关系和富集分析仍由下游应用负责。

## KEGG

独立的 KEGG 映射或 BRITE 配置发布为单表 DuckDB：

```python
from bioextract import KEGGDatabase

source = KEGGDatabase.from_brite_json("br08901.json")
source.write_duckdb("out/kegg-brite.duckdb")

published = KEGGDatabase.from_duckdb("out/kegg-brite.duckdb")
with published.connect() as connection:
    pathway_count = connection.sql("SELECT count(*) FROM pathway").fetchone()[0]
```

当多个 KEGG 产品共用一个目录时，使用最小但足以区分的限定词，例如
`kegg-mapping.duckdb` 或 `kegg-brite.duckdb`。

由化合物、反应、酶和模块组成的快照是一个多关系代谢发布产物：

```python
database = KEGGDatabase.from_metabolic_files("kegg/metabolic/2026-07")
database.write_duckdb("out/kegg.duckdb")

published = KEGGDatabase.from_duckdb("out/kegg.duckdb")
selection = published.select_ids(["CHEBI:15377"], namespace="chebi")

lf_reactions = selection.reactions()
lf_pathways = selection.pathway_memberships()
lf_unmatched = selection.unmatched_ids()

with published.connect() as connection:
    relation = connection.sql(
        """
        SELECT reaction_id, count(*) AS participant_count
        FROM reaction_participant
        GROUP BY reaction_id
        """
    )
```

领域 API 提供以反应为中心的遍历与输入血缘。`connect()` 将同一份已验证发布产物
作为由调用方管理的原生只读 DuckDB SQL 接口暴露出来。

## Reactome 与 WikiPathways

通路实体及其成员关系一同发布：

```python
from bioextract import ReactomeDatabase, WikiPathwaysDatabase

ReactomeDatabase.from_files(
    uniprot_mapping="UniProt2Reactome.txt",
    pathways="ReactomePathways.txt",
    relations="ReactomePathwaysRelation.txt",
).write_duckdb("out/reactome.duckdb")

WikiPathwaysDatabase.from_gmt(
    "wikipathways-20260510-gmt-*.gmt",
    species="Homo sapiens",
).write_duckdb("out/wikipathways.duckdb")
```

`select_ids()` 等选择方法会保留未匹配输入，并隐藏资源特有的映射连接；它们不会
计算富集统计量。WikiPathways 默认展开 glob；若要传入字面路径或路径序列，请设置
`glob=False`。构造器会先验证完整解析文件集只有一个 Collection、一个 Version，
且通路 ID 唯一，然后才应用可选的物种行过滤。

## eggNOG

提取时直接使用官方 SQLite 表示：

```python
from bioextract import EggNOGDatabase

db = EggNOGDatabase.from_sqlite(
    "eggnog.db",
    cog_functions="cog-24.fun.tab",
)
mapping = db.select_ids(["9606.ENSP00000369497"]).mappings()
```

选择操作直接查询未压缩 SQLite，无需发布衍生产物。可以传入 gzip 封装的 SQLite，
但会收到警告，并且它只会被解压到临时工作空间；若需重复使用，请先解压一次，再
传入 `.db` 文件。

## InterPro 与 Pfam

InterPro 映射，以及配置来源文件中所有可用的 Pfam 术语、交叉引用和
蛋白质-术语关系，共享一个 DuckDB 发布产物：

```python
from bioextract import InterProDatabase

db = InterProDatabase.from_mapping_files(
    protein_to_interpro="108.0/raw/protein2ipr.dat.gz",
    interpro_xml="108.0/raw/interpro.xml.gz",
)
db.write_duckdb("out/interpro.duckdb")

published = InterProDatabase.from_duckdb("out/interpro.duckdb")
with published.connect() as connection:
    print(connection.sql("SHOW TABLES").fetchall())
```

## UniProt

UniProt idmapping 仍是独立的惰性来源配置，并在 DuckDB 中发布一个 `mapping`
表：

```python
from bioextract import UniProtDatabase

UniProtDatabase.from_idmapping(
    "idmapping_selected.tab.gz",
    release_version="2026_01",
).write_duckdb(
    "out/uniprot_idmapping.duckdb",
    taxon_ids=["9606", "10090"],
)

mapping = UniProtDatabase.from_duckdb("out/uniprot_idmapping.duckdb")
human = mapping.scan_mapping(taxon_ids=["9606"])
with mapping.connect() as connection:
    print(connection.sql("SELECT count(*) FROM mapping").fetchone())
```

经过审校的 UniProtKB 是多关系 DuckDB 发布产物：

```python
UniProtDatabase.from_knowledgebase(
    entries="uniprot_sprot.dat.gz",
    canonical_sequences="uniprot_sprot.fasta.gz",
    isoform_sequences="uniprot_sprot_varsplic.fasta.gz",
    release_version="2026_01",
).write_duckdb("out/uniprot.duckdb")

db = UniProtDatabase.from_duckdb("out/uniprot.duckdb")
proteins = db.select_ids(
    ["P04637"],
    namespace="uniprot",
    taxon_ids=["9606"],
).proteins()
with db.connect() as connection:
    relation_count = connection.execute(
        "SELECT count(*) FROM protein"
    ).fetchone()[0]
```

构造器参数声明来源角色，表头和记录语法则验证其内容。路径从不提供发布版本身份。
包含所有 taxid 的 idmapping 导出必须显式设置 `allow_all_taxa=True`。

## STRING

`select_ids()` 和 `select_groups()` 封装别名解析、未匹配 ID、组间隔离和边映射：

```python
from bioextract import STRINGDatabase

selection = (
    STRINGDatabase.from_files(
        aliases="9606.protein.aliases.v12.0.txt.gz",
        links="9606.protein.links.v12.0.txt.gz",
    )
    .select_groups(
        {
            "TumorA": ["TP53", "EGFR"],
            "TumorB": ["CDK2", "TP53"],
        }
    )
    .with_min_combined_score(400)
)

lf_mapping = selection.mappings()
lf_unmapped = selection.unmatched_ids()
lf_edges = selection.edges()
```

`combined_score` 是 STRING 置信分数，而不是相互作用强度的度量。

STRING 的别名选择接受普通的非 pipe 别名文本，或完整的 UniProt
`sp|accession|entry_name` / `tr|accession|entry_name` 表示。非法的含 pipe
别名会抛出 `ValueError`；`namespace="string"` 的输入只在去除两端空白后
按 STRING protein ID 原文匹配。

## OmniPath

```python
from bioextract import OmniPathDatabase

selection = (
    OmniPathDatabase.from_files(
        enzsub="enzsub.tsv.gz",
        interactions="interactions.tsv.gz",
    )
    .select_ids(["P31749", "AKT1", "BAD"])
    .with_enzsub()
)

lf_enzsub = selection.enzsub()
lf_unmapped = selection.unmatched_ids()
```

OmniPath 的 protein 选择接受普通 protein ID 或相同的完整 UniProt pipe
表示；非法的含 pipe 调用方输入会在查找前被拒绝。

## 命名与兼容性

公共资源句柄全部使用完整的 `*Database` 名称，包括 `GODatabase`、
`ChEBIDatabase`、`RheaDatabase`、`KEGGDatabase`、`ReactomeDatabase`、
`WikiPathwaysDatabase`、`EggNOGDatabase`、`InterProDatabase`、
`UniProtDatabase`、`STRINGDatabase` 和 `OmniPathDatabase`。

从延迟加载的顶层 API 导入数据库句柄：

```python
from bioextract import ChEBIDatabase, RheaDatabase
```

对应的资源子包路径（例如 `from bioextract.rhea import RheaDatabase`）仍保持稳定。
从多个资源导入句柄时，优先使用顶层形式。

通过 `bioextract.errors` 捕获公开的运行错误类别：

```python
from bioextract.errors import CapabilityError, IntegrityError
```

数据库方法会返回或接收选择、结果、命名空间、配置和 tidy 实现类型，但这些类型
并不是稳定的包导出。不要依赖其深层模块路径来获得兼容性。

本库不提供缩写的 `*Db` 别名、旧版分数过滤器名称或目录写入器。请直接使用
`with_min_combined_score()` 和 `write_duckdb()`。

表名、视图名和生成列名使用单数形式的 `snake_case`。除非为了使官方二维数据可
查询而必须进行最小的确定性映射，否则保留原始官方表头。任何此类映射都会记录
在嵌入式来源信息中。

带版本的 CephFS 约定为 `tidy/data.duckdb`。调用方可以使用其他文件名；文件名
绝不是 schema 身份或兼容性标识。机器可读的身份来自嵌入式元数据。

## 开发

- 文档索引位于 [docs/README.md](docs/README.md)。
- 测试分层和 fixture 所有权定义在
  [docs/testing/README.md](docs/testing/README.md)。
- `pdm run format`
- `pdm run lint`
- `pdm run typecheck`
- `pdm run test-unit`
- `pdm run test-contract`
- `pdm run test-integration`
- `pdm run test`
- `pdm run test-smoke` 只运行显式配置的主机发布产物测试。
- `pdm run precommit` 会应用格式化和 lint 修复，然后运行严格类型检查与完整的
  自包含测试套件。

自包含测试默认将 DuckDB、Polars 和 Rayon 驱动的任务限制为四个线程。在受限
主机上共享资源时，请设置 `BIOEXTRACT_TEST_THREADS=1`。

构建发布产物时，请在导入 Polars 或 bioextract 之前设置
`POLARS_MAX_THREADS`。它同时限制 Polars 执行和 bioextract 管理的 DuckDB 发布
连接。

## 发布

- `.github/workflows/py-ci.yml` 运行测试与构建检查。
- `.github/workflows/publish.yml` 发布规范的 PEP 440 标签。
- `pypi` 环境应配置 PyPI 可信发布。
