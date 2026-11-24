"""Unit harmonization: everything to nM, then to p-scale (pIC50, pKi, ...).

ChEMBL already provides ``standard_value``/``standard_units`` (usually nM)
and ``pchembl_value``, but students still hit mixed units (``uM``, ``mM``,
``ug.mL-1``, ``%``), values without units, and non-numeric values. These
helpers convert what *can* be converted safely and flag the rest instead of
silently mixing scales.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

# Multiplicative factor to convert a concentration in <unit> to nM.
MOLAR_TO_NM: dict[str, float] = {
    "fM": 1e-6,
    "pM": 1e-3,
    "nM": 1.0,
    "uM": 1e3,
    "µM": 1e3,
    "μM": 1e3,
    "mM": 1e6,
    "M": 1e9,
    "nmol/L": 1.0,
    "nmol.L-1": 1.0,
    "umol/L": 1e3,
    "umol.L-1": 1e3,
    "mmol/L": 1e6,
    "mol/L": 1e9,
    "10'-9M": 1.0,
    "10'-6M": 1e3,
}

# Activity types whose value is a concentration that can be put on a -log10 molar scale.
P_SCALE_TYPES = frozenset(
    {"IC50", "EC50", "Ki", "Kd", "AC50", "GI50", "XC50", "ED50", "Potency", "IC90", "MIC"}
)


def to_nM(value: float | str | None, units: str | None) -> float:
    """Convert a concentration to nM; returns ``nan`` when it is not convertible."""
    if value is None or units is None:
        return math.nan
    try:
        v = float(value)
    except (TypeError, ValueError):
        return math.nan
    factor = MOLAR_TO_NM.get(str(units).strip())
    if factor is None or not math.isfinite(v):
        return math.nan
    return v * factor


def nM_to_p(value_nM: float) -> float:
    """-log10(molar) from nM, e.g. 10 nM -> 8.0. Non-positive values give ``nan``."""
    try:
        v = float(value_nM)
    except (TypeError, ValueError):
        return math.nan
    if not math.isfinite(v) or v <= 0:
        return math.nan
    return 9.0 - math.log10(v)


def p_to_nM(p_value: float) -> float:
    """Inverse of :func:`nM_to_p`."""
    return 10 ** (9.0 - float(p_value))


def harmonize_units(
    df: pd.DataFrame,
    value_col: str = "standard_value",
    units_col: str = "standard_units",
    type_col: str = "standard_type",
) -> pd.DataFrame:
    """Add ``value_nM``, ``p_activity`` and ``unit_status`` columns.

    ``p_activity`` uses ChEMBL's ``pchembl_value`` when present (it is curated
    and only set for clean ``=`` measurements), otherwise -log10(M) computed
    from ``value_nM``. ``unit_status`` is ``"ok"``, ``"unconvertible_unit"``,
    ``"missing_value"`` or ``"non_positive"``.
    """
    out = df.copy()
    if out.empty:
        for col in ("value_nM", "p_activity", "unit_status"):
            out[col] = pd.Series(dtype=float if col != "unit_status" else object)
        return out
    values = pd.to_numeric(out[value_col], errors="coerce")
    factors = out[units_col].map(lambda u: MOLAR_TO_NM.get(str(u).strip()) if pd.notna(u) else None)
    factors = pd.to_numeric(factors, errors="coerce")
    out["value_nM"] = values * factors

    computed = out["value_nM"].map(nM_to_p)
    if type_col in out:
        computed = computed.where(out[type_col].isin(P_SCALE_TYPES))
    if "pchembl_value" in out:
        pchembl = pd.to_numeric(out["pchembl_value"], errors="coerce")
        out["p_activity"] = pchembl.fillna(computed)
    else:
        out["p_activity"] = computed

    status = np.where(values.isna(), "missing_value",
             np.where(factors.isna(), "unconvertible_unit",
             np.where(values <= 0, "non_positive", "ok")))
    out["unit_status"] = status
    return out
