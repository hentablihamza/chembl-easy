# Offline mode: the ChEMBL SQLite dump

The REST API is ideal for one target or a few hundred compounds. For
**whole-database work** (all kinases, every IC50 in ChEMBL, repeated large
pulls on a cluster without internet) download the SQLite release once and
query it locally.

| | REST API (`chembl-easy`) | SQLite dump |
| --- | --- | --- |
| Setup | none | ~5.8 GB download (`chembl_37_sqlite.tar.gz`), more once extracted |
| One target (e.g. EGFR) | minutes the first time, seconds from cache | seconds |
| All of ChEMBL | impractical | the right tool |
| Always latest release | yes | until the next release |

## 1. Download

```bash
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/checksums.txt
sha256sum -c checksums.txt --ignore-missing
tar -xzf chembl_37_sqlite.tar.gz      # -> chembl_37/chembl_37_sqlite/chembl_37.db
```

Replace `37` with the current release number (`chembl-easy status` prints it).

## 2. Same columns as the API, so the rest of chembl-easy still works

The query below returns the same column names that `fetch_activities`
returns, so `harmonize_units` and `aggregate_by_compound` work unchanged.
Table and column names follow the ChEMBL 37 schema documentation
(`schema_documentation.txt` on the FTP site).

```python
import sqlite3
import pandas as pd
from chembl_easy import harmonize_units, aggregate_by_compound

SQL = """
SELECT act.activity_id,
       md.chembl_id            AS molecule_chembl_id,
       pmd.chembl_id           AS parent_molecule_chembl_id,
       md.pref_name            AS molecule_pref_name,
       cs.canonical_smiles,
       act.standard_type, act.standard_relation, act.standard_value, act.standard_units,
       act.pchembl_value, act.activity_comment, act.data_validity_comment, act.potential_duplicate,
       a.chembl_id             AS assay_chembl_id,
       a.assay_type, a.description AS assay_description, a.confidence_score, a.variant_id,
       td.chembl_id            AS target_chembl_id,
       td.organism             AS target_organism,
       d.chembl_id             AS document_chembl_id,
       d.year                  AS document_year
FROM activities act
JOIN assays a              ON act.assay_id = a.assay_id
JOIN target_dictionary td  ON a.tid = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh  ON act.molregno = mh.molregno
LEFT JOIN molecule_dictionary pmd ON mh.parent_molregno = pmd.molregno
LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
LEFT JOIN docs d           ON act.doc_id = d.doc_id
WHERE td.chembl_id = ?
  AND act.standard_type = 'IC50'
  AND act.standard_relation = '='
  AND a.assay_type = 'B'
  AND a.confidence_score >= 9
  AND a.variant_id IS NULL                 -- wild type only
  AND act.data_validity_comment IS NULL    -- not flagged by curators
  AND COALESCE(act.potential_duplicate, 0) = 0
  AND cs.canonical_smiles IS NOT NULL
"""

with sqlite3.connect("chembl_37/chembl_37_sqlite/chembl_37.db") as con:
    raw = pd.read_sql_query(SQL, con, params=("CHEMBL203",))

clean = aggregate_by_compound(harmonize_units(raw).dropna(subset=["p_activity"]))
```

Note the one difference: the API exposes `assay_variant_mutation` (text),
while the database stores `assays.variant_id` (a foreign key to
`variant_sequences`); the query above filters on it directly.

Results can differ slightly from the API if the two are on different ChEMBL
releases, so always report the release you used.
