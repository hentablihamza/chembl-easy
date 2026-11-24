import pytest

from chembl_easy.chem import InvalidSmilesError
from chembl_easy.molecules import flatten_molecule, get_molecule, search_molecules, similarity_search, substructure_search

from conftest import FakeResponse, paged


def mol(cid, name=None, smiles="CCO", sim=None):
    rec = {
        "molecule_chembl_id": cid,
        "pref_name": name,
        "max_phase": "4.0" if name else None,
        "molecule_type": "Small molecule",
        "molecule_structures": {"canonical_smiles": smiles, "standard_inchi_key": "KEY", "molfile": "..."},
        "molecule_properties": {"full_mwt": "46.07", "alogp": "-0.03", "hba": 1, "hbd": 1, "psa": "20.23",
                                "rtb": 0, "num_ro5_violations": 0, "qed_weighted": "0.41"},
    }
    if sim is not None:
        rec["similarity"] = sim
    return rec


IMATINIB = mol("CHEMBL941", "IMATINIB", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1")


def test_flatten_molecule_drops_molfile():
    row = flatten_molecule(mol("CHEMBL545", "ETHANOL", sim="88.5"))
    assert row["canonical_smiles"] == "CCO"
    assert row["full_mwt"] == "46.07"
    assert row["similarity"] == 88.5
    assert "molfile" not in row
    assert flatten_molecule({"molecule_chembl_id": "X", "molecule_structures": None})["canonical_smiles"] is None


def test_get_molecule(make_client):
    client, session = make_client(lambda ep, p: IMATINIB if ep == "molecule/CHEMBL941.json" else FakeResponse(404))
    assert get_molecule("chembl941", client=client)["pref_name"] == "IMATINIB"


def test_search_by_trade_name_uses_synonyms(make_client):
    def handler(ep, p):
        if "molecule_synonyms__molecule_synonym__iexact" in p:
            assert p["molecule_synonyms__molecule_synonym__iexact"] == "Gleevec"
            return paged("molecules", [IMATINIB], p)
        if "pref_name__iexact" in p:
            return paged("molecules", [], p)
        raise AssertionError("full-text search should not be needed")

    client, _ = make_client(handler)
    df = search_molecules("Gleevec", client=client)
    assert list(df["molecule_chembl_id"]) == ["CHEMBL941"]
    assert df.loc[0, "full_mwt"] == pytest.approx(46.07)


def test_search_dedupes_and_falls_back_to_full_text(make_client):
    def handler(ep, p):
        if ep == "molecule/search.json":
            return paged("molecules", [IMATINIB, IMATINIB, mol("CHEMBL1642", "IMATINIB MESYLATE")], p)
        return paged("molecules", [], p)

    client, _ = make_client(handler)
    df = search_molecules("imatinib kinase", client=client)
    assert list(df["molecule_chembl_id"]) == ["CHEMBL941", "CHEMBL1642"]


def test_search_by_chembl_id(make_client):
    client, _ = make_client(lambda ep, p: IMATINIB if ep == "molecule/CHEMBL941.json" else FakeResponse(404))
    assert len(search_molecules("CHEMBL941", client=client)) == 1
    assert search_molecules("CHEMBL0", client=client).empty


def test_similarity_sends_smiles_as_parameter(make_client):
    erlotinib = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

    def handler(ep, p):
        assert ep == "similarity.json"  # never in the URL path, so '#' is safe
        assert p["smiles"] == erlotinib and p["similarity"] == 80
        return paged("molecules", [mol("A", sim="90.1"), mol("CHEMBL553", "ERLOTINIB", sim="100")], p)

    client, session = make_client(handler)
    df = similarity_search(erlotinib, threshold=80, client=client)
    assert list(df["molecule_chembl_id"]) == ["CHEMBL553", "A"]
    assert "#" not in session.calls[0][0]


def test_similarity_validation(make_client):
    client, session = make_client(lambda ep, p: AssertionError("should not be called"))
    with pytest.raises(ValueError, match="40 and 100"):
        similarity_search("CCO", threshold=30, client=client)
    with pytest.raises(InvalidSmilesError):
        similarity_search("C1CC(", client=client)
    assert session.calls == []


def test_substructure(make_client):
    records = [mol(f"CHEMBL{i}") for i in range(30)]
    client, _ = make_client(lambda ep, p: paged("molecules", records, p))
    df = substructure_search("c1ccc2ncncc2c1", client=client, max_records=12)
    assert len(df) == 12
    with pytest.raises(InvalidSmilesError):
        substructure_search("((", client=client)
