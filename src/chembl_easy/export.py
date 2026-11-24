"""Write tables to CSV, TSV or Parquet, chosen from the file extension."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def save_table(df: pd.DataFrame, path: str | Path) -> Path:
    """Save ``df`` as ``.csv``, ``.tsv``/``.txt`` or ``.parquet``; parent dirs are created."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        df.to_csv(path, index=False)
    elif suffix in {".tsv", ".txt"}:
        df.to_csv(path, index=False, sep="\t")
    elif suffix in {".parquet", ".pq"}:
        try:
            df.to_parquet(path, index=False)
        except ImportError as exc:
            raise ImportError(
                "Parquet export needs pyarrow: pip install 'chembl-easy[parquet]'"
            ) from exc
    else:
        raise ValueError(f"Unsupported extension {suffix!r}; use .csv, .tsv or .parquet")
    return path


def ml_ready(df: pd.DataFrame, threshold: float | None = None) -> pd.DataFrame:
    """Minimal ``compound_id, smiles, p_activity[, active]`` table for modelling."""
    out = df[["compound_id", "smiles", "p_activity"]].copy()
    if threshold is not None:
        out["active"] = (out["p_activity"] >= threshold).astype(int)
    return out.reset_index(drop=True)
