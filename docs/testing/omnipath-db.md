# OmniPath Test Standard

OmniPath tests use temporary local TSV fixtures and cover the direct-adapter
contract without network or CephFS access.

- source scans preserve official headers and order;
- optional `enzsub` and `interactions` resources fail clearly when absent;
- `protein` namespace validation is explicit;
- single and grouped selections preserve evidence fields and unmatched IDs;
- grouped extraction resolves shared identifiers once and retains empty groups;
- plain and gzip inputs follow the same relation semantics; and
- malformed headers are rejected at the source-profile boundary.

Run the focused suite with:

```console
PYTHONPATH=src python -m pytest tests/integration/omnipath/test_database.py -q
```

Real OmniPath snapshots are an explicit external smoke concern and are not
required for the hermetic repository gate.
