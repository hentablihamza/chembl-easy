# Example data

`egfr/` was produced by `python examples/egfr_example.py` on 24 November 2026
against the live ChEMBL REST API, release **ChEMBL_37**.

| File | Content |
| --- | --- |
| `egfr_ic50_compounds.csv` | 6,871 parent compounds with median pIC50 and per-compound statistics |
| `egfr_ic50_ml.csv` | `compound_id, smiles, p_activity, active` (active = pIC50 >= 7) |
| `summary.json` | counts, pIC50 distribution, raw units, filtering funnel, request count |
| `funnel.md` | filtering funnel as a Markdown table |
| `p_activity_hist.png` | pIC50 histogram |

Filters: target CHEMBL203 (EGFR, *Homo sapiens*, single protein), standard type
IC50, relation `=`, binding assays (B), assay confidence score 9, wild-type
assays only, records without `data_validity_comment`, no `potential_duplicate`.

ChEMBL data is provided by EMBL-EBI under Creative Commons
Attribution-ShareAlike 3.0 Unported; ChEMBL IDs are preserved. Cite:
Mendez D. *et al.*, *Nucleic Acids Res.* 2019, 47(D1), D930–D940,
doi:10.1093/nar/gky1075.
