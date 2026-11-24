"""Collapse repeated measurements into one ML-ready row per compound.

The same compound is often measured many times (different papers, assays or
salt forms). Training on raw rows over-weights popular compounds and leaks
duplicates across train/test splits. :func:`aggregate_by_compound` groups by
parent compound, summarizes the p-activity, and flags compounds whose
measurements disagree by more than ``max_range`` log units.
"""

from __future__ import annotations

import pandas as pd

from chembl_easy.chem import parent_smiles

AGG_COLUMNS = [
    "compound_id",
    "smiles",
    "p_activity",
    "p_activity_mean",
    "p_activity_std",
    "p_activity_min",
    "p_activity_max",
    "n_measurements",
    "n_assays",
    "n_documents",
    "standard_types",
    "first_year",
    "high_variance",
]


def aggregate_by_compound(
    df: pd.DataFrame,
    value_col: str = "p_activity",
    method: str = "median",
    max_range: float = 1.0,
    drop_high_variance: bool = False,
    strip_salts: bool = True,
) -> pd.DataFrame:
    """One row per parent compound with a consensus ``p_activity``.

    Compounds are grouped by ``parent_molecule_chembl_id`` (so salt forms
    merge), falling back to ``molecule_chembl_id``. ``smiles`` is the
    salt-stripped, neutralized RDKit canonical SMILES when ``strip_salts``.
    ``high_variance`` is True when max - min of the measurements exceeds
    ``max_range`` log units.
    """
    if method not in {"median", "mean"}:
        raise ValueError("method must be 'median' or 'mean'")
    data = df[pd.to_numeric(df[value_col], errors="coerce").notna()].copy()
    data[value_col] = pd.to_numeric(data[value_col])
    if data.empty:
        return pd.DataFrame(columns=AGG_COLUMNS)

    if "parent_molecule_chembl_id" in data:
        data["compound_id"] = data["parent_molecule_chembl_id"].fillna(data["molecule_chembl_id"])
    else:
        data["compound_id"] = data["molecule_chembl_id"]
    if strip_salts:
        data["smiles"] = data["canonical_smiles"].map(parent_smiles)
    else:
        data["smiles"] = data["canonical_smiles"]
    data = data[data["smiles"].notna()]

    def _join(values: pd.Series) -> str:
        return ";".join(sorted({str(v) for v in values.dropna()}))

    g = data.groupby("compound_id", sort=True)
    out = pd.DataFrame(
        {
            "smiles": g["smiles"].agg(lambda s: s.mode().iloc[0]),
            "p_activity": g[value_col].median() if method == "median" else g[value_col].mean(),
            "p_activity_mean": g[value_col].mean(),
            "p_activity_std": g[value_col].std(ddof=1).fillna(0.0),
            "p_activity_min": g[value_col].min(),
            "p_activity_max": g[value_col].max(),
            "n_measurements": g.size(),
            "n_assays": g["assay_chembl_id"].nunique() if "assay_chembl_id" in data else 1,
            "n_documents": g["document_chembl_id"].nunique() if "document_chembl_id" in data else 1,
            "standard_types": g["standard_type"].agg(_join) if "standard_type" in data else "",
            "first_year": pd.to_numeric(data["document_year"], errors="coerce").groupby(data["compound_id"]).min()
            if "document_year" in data
            else pd.NA,
        }
    )
    out["high_variance"] = (out["p_activity_max"] - out["p_activity_min"]) > max_range
    if drop_high_variance:
        out = out[~out["high_variance"]]
    out = out.reset_index()
    for col in ("p_activity", "p_activity_mean", "p_activity_std", "p_activity_min", "p_activity_max"):
        out[col] = out[col].round(3)
    return out[AGG_COLUMNS]


def add_activity_class(df: pd.DataFrame, threshold: float = 6.0, col: str = "p_activity") -> pd.DataFrame:
    """Add a binary ``active`` column (p_activity >= threshold, default 1 uM)."""
    out = df.copy()
    out["active"] = (out[col] >= threshold).astype(int)
    return out
