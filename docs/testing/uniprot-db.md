# UniProtDatabase Test Standard

The focused suite verifies:

- raw plain/gzip, Parquet, and legacy hive idmapping scans;
- exact 22-column schema validation and scoped eager reads;
- atomic one-table idmapping DuckDB publication with no sidecar;
- metadata-v1 exactness and idmapping/knowledgebase profile discrimination;
- source/reopened mapping parity, taxon scoping, all-taxa safety, fresh
  caller-owned read-only connections, and native SQL;
- role declaration independent of basename and compression suffix;
- strict reviewed DAT grammar and mandatory `ID`, `AC`, `OX`, and `SQ` facts;
- SQ length, molecular weight, CRC64, and sequence parity;
- exact canonical FASTA equality and varsplic-to-DAT isoform resolution;
- multi-IsoId products retain one main product row plus ordered old-ID rows;
- repeated External/owned product contexts remain distinct, while varsplic
  materializes only the unique Alternative owner;
- metadata v1 resource/source profile fields and optional release identity;
- atomic DuckDB publication, inventory parity, and read-only native SQL;
- primary/secondary/isoform-ID namespace selection, grouping, core extractors,
  and unmatched IDs.

The canonical repository gate is `pdm run check`. The 2026_01 full snapshot
smoke is non-publishing and must not write CephFS.
