import json

import pandas as pd
import pytest

from chembl_easy import cli
from chembl_easy.activities import ActivityQuery
from chembl_easy.client import ChEMBLClient
from chembl_easy.export import ml_ready, save_table
from chembl_easy.pipeline import build_target_dataset

from conftest import FakeResponse, FakeSession, activity, paged, target

ACTS = [
    activity(1, canonical_smiles="Oc1ccccc1", pchembl_value="7.0"),
    activity(2, canonical_smiles="Oc1ccccc1", pchembl_value="8.0", assay_chembl_id="A2"),
    activity(3, molecule_chembl_id="CHEMBL3s", parent_molecule_chembl_id="CHEMBL3",
             canonical_smiles="CCN.Cl", pchembl_value=None, standard_value="1", standard_units="uM"),
    activity(4, canonical_smiles="CCCl", pchembl_value=None, standard_value="5", standard_units="ug.mL-1"),
]


def handler(ep, p):
    if ep == "target.json":
        return paged("targets", [target("CHEMBL203", "Epidermal growth factor receptor", "SINGLE PROTEIN")], p)
    if ep == "activity.json":
        return paged("activities", ACTS, p)
    if ep == "status.json":
        return {"chembl_db_version": "ChEMBL_37"}
    if ep == "molecule.json" or ep == "molecule/search.json":
        return paged("molecules", [{"molecule_chembl_id": "CHEMBL941", "pref_name": "IMATINIB",
                                    "molecule_structures": {"canonical_smiles": "CCO"}}], p)
    if ep == "similarity.json":
        return paged("molecules", [{"molecule_chembl_id": "CHEMBL545", "similarity": "100",
                                    "molecule_structures": {"canonical_smiles": "CCO"}}], p)
    if ep == "substructure.json":
        return paged("molecules", [{"molecule_chembl_id": "CHEMBL545"}], p)
    return FakeResponse(404)


@pytest.fixture
def fake_cli(monkeypatch):
    session = FakeSession(handler)

    def _client(args):
        return ChEMBLClient(cache_dir=None, session=session, sleep=lambda s: None, min_interval=0)

    monkeypatch.setattr(cli, "_client", _client)
    return session


def test_save_table_formats(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    for name in ["t.csv", "sub/t.tsv", "t.parquet"]:
        path = save_table(df, tmp_path / name)
        back = pd.read_parquet(path) if name.endswith("parquet") else pd.read_csv(path, sep="\t" if name.endswith("tsv") else ",")
        pd.testing.assert_frame_equal(back, df)
    with pytest.raises(ValueError):
        save_table(df, tmp_path / "t.xlsx")


def test_ml_ready():
    df = pd.DataFrame({"compound_id": ["a", "b"], "smiles": ["C", "CC"], "p_activity": [6.5, 7.5], "x": [1, 2]})
    assert list(ml_ready(df).columns) == ["compound_id", "smiles", "p_activity"]
    assert list(ml_ready(df, threshold=7)["active"]) == [0, 1]


def test_build_target_dataset(make_client):
    client, _ = make_client(handler)
    ds = build_target_dataset("EGFR", query=ActivityQuery(), client=client)
    assert ds.target["target_chembl_id"] == "CHEMBL203"
    assert ds.chembl_release == "ChEMBL_37"
    clean = ds.clean.set_index("compound_id")
    assert set(clean.index) == {"CHEMBL1", "CHEMBL2", "CHEMBL3"}  # ug/mL record is not convertible
    assert clean.loc["CHEMBL3", "smiles"] == "CCN"
    assert clean.loc["CHEMBL3", "p_activity"] == pytest.approx(6.0)
    steps = dict(ds.report.steps)
    assert steps["Convertible to p-scale (pIC50/pKi/...)"] == 3
    assert steps["Unique parent compounds"] == 3
    s = ds.summary()
    assert s["compounds"] == 3 and s["actives_p_ge_7"] == 2


def test_cli_activities_writes_files(fake_cli, tmp_path, capsys):
    out, raw, summ = tmp_path / "c.csv", tmp_path / "raw.parquet", tmp_path / "s.json"
    rc = cli.main(["activities", "EGFR", "-o", str(out), "--raw-output", str(raw), "--summary", str(summ),
                   "--ml-ready", "--active-threshold", "7", "-q"])
    assert rc == 0
    df = pd.read_csv(out)
    assert list(df.columns) == ["compound_id", "smiles", "p_activity", "active"]
    assert len(pd.read_parquet(raw)) == 4
    summary = json.loads(summ.read_text())
    assert summary["target_chembl_id"] == "CHEMBL203" and summary["compounds"] == 3
    assert "CHEMBL203" in capsys.readouterr().out
    params = [p for url, p in fake_cli.calls if url.endswith("activity.json")][0]
    assert params["standard_type__in"] == "IC50" and params["standard_relation__in"] == "="


def test_cli_activities_options(fake_cli, capsys):
    rc = cli.main(["activities", "EGFR", "--types", "IC50", "Ki", "--any-relation", "--assay-types", "B", "F",
                   "--pchembl-only", "-q"])
    assert rc == 0
    params = [p for url, p in fake_cli.calls if url.endswith("activity.json")][0]
    assert params["standard_type__in"] == "IC50,Ki"
    assert "standard_relation__in" not in params
    assert params["assay_type__in"] == "B,F" and params["pchembl_value__isnull"] == "false"


def test_cli_target_molecule_similar_substructure(fake_cli, capsys):
    assert cli.main(["target", "EGFR"]) == 0
    assert "CHEMBL203" in capsys.readouterr().out
    assert cli.main(["molecule", "imatinib"]) == 0
    assert "CHEMBL941" in capsys.readouterr().out
    assert cli.main(["similar", "CCO", "--threshold", "90"]) == 0
    assert "CHEMBL545" in capsys.readouterr().out
    assert cli.main(["substructure", "CCO", "--max-records", "1"]) == 0


def test_cli_errors(fake_cli, capsys):
    assert cli.main(["similar", "C1CC("]) == 2
    assert "cannot parse" in capsys.readouterr().err
    assert cli.main(["similar", "CCO", "--threshold", "20"]) == 2


def test_cli_status_and_clear_cache(fake_cli, tmp_path, capsys):
    assert cli.main(["status"]) == 0
    assert "ChEMBL_37" in capsys.readouterr().out
    assert cli.main(["--cache-dir", str(tmp_path), "clear-cache"]) == 0


def test_cli_version(capsys):
    with pytest.raises(SystemExit):
        cli.main(["--version"])
    assert "chembl-easy" in capsys.readouterr().out
