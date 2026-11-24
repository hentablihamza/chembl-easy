import math

import pandas as pd
import pytest

from chembl_easy.units import harmonize_units, nM_to_p, p_to_nM, to_nM


@pytest.mark.parametrize(
    "value,units,expected",
    [(100, "nM", 100), (1, "uM", 1000), (1, "µM", 1000), (2, "mM", 2e6), (5, "pM", 0.005), ("3", "nmol/L", 3)],
)
def test_to_nM(value, units, expected):
    assert to_nM(value, units) == pytest.approx(expected)


@pytest.mark.parametrize("value,units", [(1, "ug.mL-1"), (1, "%"), (None, "nM"), ("abc", "nM"), (1, None)])
def test_to_nM_unconvertible(value, units):
    assert math.isnan(to_nM(value, units))


def test_p_scale():
    assert nM_to_p(10) == pytest.approx(8.0)
    assert nM_to_p(1000) == pytest.approx(6.0)
    assert math.isnan(nM_to_p(0))
    assert math.isnan(nM_to_p(-5))
    assert p_to_nM(7) == pytest.approx(100)


def test_harmonize_prefers_pchembl_and_flags_units():
    df = pd.DataFrame(
        {
            "standard_type": ["IC50", "IC50", "IC50", "Ki", "Inhibition"],
            "standard_value": ["100", "1", "5", "0", "50"],
            "standard_units": ["nM", "uM", "ug.mL-1", "nM", "%"],
            "pchembl_value": ["7.10", None, None, None, None],
        }
    )
    out = harmonize_units(df)
    assert out.loc[0, "p_activity"] == pytest.approx(7.10)  # curated pChEMBL wins
    assert out.loc[1, "value_nM"] == pytest.approx(1000)
    assert out.loc[1, "p_activity"] == pytest.approx(6.0)
    assert math.isnan(out.loc[2, "p_activity"])
    assert list(out["unit_status"]) == ["ok", "ok", "unconvertible_unit", "non_positive", "unconvertible_unit"]
    assert math.isnan(out.loc[4, "p_activity"])


def test_harmonize_empty():
    out = harmonize_units(pd.DataFrame(columns=["standard_value", "standard_units", "standard_type"]))
    assert {"value_nM", "p_activity", "unit_status"} <= set(out.columns)
