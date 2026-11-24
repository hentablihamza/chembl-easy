import pandas as pd
import pytest

from chembl_easy.aggregate import add_activity_class, aggregate_by_compound
from chembl_easy.chem import InvalidSmilesError, inchikey, parent_smiles, validate_smiles


def test_validate_smiles():
    assert validate_smiles("OC1=CC=CC=C1") == "Oc1ccccc1"
    for bad in ["C1CC(", "", "   ", "XYZ"]:
        with pytest.raises(InvalidSmilesError):
            validate_smiles(bad)


def test_parent_smiles_strips_salts_and_neutralizes():
    assert parent_smiles("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]") == "CC(=O)Oc1ccccc1C(=O)O"
    assert parent_smiles("CS(=O)(=O)O.Cc1ccccc1N") == "Cc1ccccc1N"
    assert parent_smiles("not a smiles") is None
    assert parent_smiles(None) is None


def test_inchikey():
    assert inchikey("CC(=O)Oc1ccccc1C(=O)O") == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert inchikey("garbage") is None


def _rows():
    return pd.DataFrame(
        {
            "molecule_chembl_id": ["M1", "M1_salt", "M1", "M2", "M3", "M4"],
            "parent_molecule_chembl_id": ["M1", "M1", "M1", "M2", "M3", None],
            "canonical_smiles": ["Oc1ccccc1", "Oc1ccccc1.Cl", "Oc1ccccc1", "CCO", "CCN", "CCC"],
            "p_activity": [6.0, 6.4, 8.0, 5.0, None, 7.0],
            "assay_chembl_id": ["A1", "A2", "A3", "A1", "A1", "A1"],
            "document_chembl_id": ["D1", "D1", "D2", "D1", "D1", "D1"],
            "standard_type": ["IC50", "IC50", "Ki", "IC50", "IC50", "IC50"],
            "document_year": [2001, 2005, 1999, 2010, 2010, 2012],
        }
    )


def test_aggregate_by_parent_compound():
    out = aggregate_by_compound(_rows()).set_index("compound_id")
    assert set(out.index) == {"M1", "M2", "M4"}  # M3 has no value; M4 falls back to molecule id
    m1 = out.loc["M1"]
    assert m1["n_measurements"] == 3
    assert m1["p_activity"] == pytest.approx(6.4)  # median of 6.0, 6.4, 8.0
    assert m1["p_activity_mean"] == pytest.approx(6.8)
    assert m1["n_assays"] == 3 and m1["n_documents"] == 2
    assert m1["standard_types"] == "IC50;Ki"
    assert m1["first_year"] == 1999
    assert m1["smiles"] == "Oc1ccccc1"  # salt stripped
    assert bool(m1["high_variance"])  # range 2.0 > 1.0
    assert not bool(out.loc["M2", "high_variance"])
    assert out.loc["M2", "p_activity_std"] == 0.0


def test_aggregate_mean_and_drop_high_variance():
    out = aggregate_by_compound(_rows(), method="mean", drop_high_variance=True)
    assert "M1" not in set(out["compound_id"])
    with pytest.raises(ValueError):
        aggregate_by_compound(_rows(), method="mode")


def test_aggregate_without_salt_stripping():
    out = aggregate_by_compound(_rows(), strip_salts=False).set_index("compound_id")
    assert out.loc["M1", "smiles"] == "Oc1ccccc1"  # most common raw SMILES


def test_aggregate_empty():
    df = _rows().assign(p_activity=None)
    out = aggregate_by_compound(df)
    assert out.empty and "p_activity" in out.columns


def test_add_activity_class():
    df = pd.DataFrame({"p_activity": [5.9, 6.0, 7.5]})
    assert list(add_activity_class(df)["active"]) == [0, 1, 1]
    assert list(add_activity_class(df, threshold=7)["active"]) == [0, 0, 1]
