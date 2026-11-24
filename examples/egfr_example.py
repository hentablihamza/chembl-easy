"""Worked example: a clean EGFR (CHEMBL203) IC50 dataset from a gene symbol.

Run:  python examples/egfr_example.py
Writes examples/egfr/{egfr_ic50_compounds.csv, egfr_ic50_ml.csv, summary.json,
funnel.md, p_activity_hist.png}. The first run downloads ~18k activity records
(a few minutes); later runs are served from the on-disk cache.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from chembl_easy import ActivityQuery, ChEMBLClient, build_target_dataset, find_targets, ml_ready, save_table

OUT = Path(__file__).parent / "egfr"
logging.basicConfig(level=logging.INFO, format="%(message)s")


def main() -> None:
    client = ChEMBLClient()

    # 1. "EGFR" is ambiguous - look at every candidate before downloading.
    candidates = find_targets("EGFR", client=client)
    print(candidates[["target_chembl_id", "pref_name", "target_type", "organism"]].head(8).to_string(index=False))
    print(f"{len(candidates)} ChEMBL targets match 'EGFR'\n")

    # 2. Binding IC50s, exact '=' values, assay confidence 9, wild type, no flagged/duplicate records.
    query = ActivityQuery(standard_types=["IC50"], relations=["="], assay_types=["B"], min_confidence=9)
    ds = build_target_dataset("EGFR", query=query, client=client)

    raw_units = ds.raw["standard_units"].fillna("(none)").value_counts().to_dict()
    summary = {**ds.summary(), "standard_units_in_raw": raw_units,
               "funnel": ds.report.steps, "http_requests": client.n_requests,
               "cache_hits": client.n_cache_hits}

    OUT.mkdir(exist_ok=True)
    save_table(ds.clean, OUT / "egfr_ic50_compounds.csv")
    save_table(ml_ready(ds.clean, threshold=7.0), OUT / "egfr_ic50_ml.csv")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "funnel.md").write_text(ds.report.to_markdown() + "\n", encoding="utf-8")
    print(ds.report.to_markdown())
    print(json.dumps(summary, indent=2))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.hist(ds.clean["p_activity"], bins=40, color="#3b6ea5", edgecolor="white")
    ax.axvline(7.0, color="#c0392b", ls="--", lw=1, label="pIC50 = 7 (100 nM)")
    ax.set_xlabel("pIC50 (median per compound)")
    ax.set_ylabel("compounds")
    ax.set_title(f"EGFR (CHEMBL203), {len(ds.clean):,} compounds, {ds.chembl_release}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "p_activity_hist.png", dpi=150)


if __name__ == "__main__":
    main()
