import pytest

from chembl_easy.targets import classify_query, find_targets, resolve_target

from conftest import paged, target

EGFR_TARGETS = [
    target("CHEMBL4523747", "EGFR/PPP1CA", "PROTEIN-PROTEIN INTERACTION", accessions=("P00533", "P62136")),
    target("CHEMBL2363049", "Epidermal growth factor receptor", "PROTEIN FAMILY", accessions=("P00533", "P04626")),
    target("CHEMBL3608", "Epidermal growth factor receptor", "SINGLE PROTEIN", "Mus musculus", ("Q01279",)),
    target("CHEMBL203", "Epidermal growth factor receptor", "SINGLE PROTEIN"),
]


@pytest.mark.parametrize(
    "query,kind",
    [("CHEMBL203", "chembl_id"), ("chembl203", "chembl_id"), ("P00533", "uniprot"), ("A0A024R161", "uniprot"),
     ("EGFR", "text"), ("beta-secretase 1", "text"), ("p00533", "text")],
)
def test_classify_query(query, kind):
    assert classify_query(query) == kind


def test_gene_symbol_ranks_human_single_protein_first(make_client):
    def handler(ep, p):
        assert ep == "target.json"
        assert p["target_components__target_component_synonyms__component_synonym__iexact"] == "EGFR"
        return paged("targets", EGFR_TARGETS, p)

    client, _ = make_client(handler)
    df = find_targets("EGFR", client=client)
    assert list(df["target_chembl_id"]) == ["CHEMBL203", "CHEMBL2363049", "CHEMBL4523747", "CHEMBL3608"]
    assert df.loc[0, "accessions"] == "P00533"
    assert set(df["match"]) == {"gene_symbol"}


def test_organism_preference_changes_ranking(make_client):
    client, _ = make_client(lambda ep, p: paged("targets", EGFR_TARGETS, p))
    df = find_targets("EGFR", client=client, organism="Mus musculus")
    assert df.loc[0, "target_chembl_id"] == "CHEMBL3608"


def test_single_protein_only(make_client):
    client, _ = make_client(lambda ep, p: paged("targets", EGFR_TARGETS, p))
    df = find_targets("EGFR", client=client, single_protein_only=True)
    assert set(df["target_type"]) == {"SINGLE PROTEIN"}


def test_free_text_fallback(make_client):
    def handler(ep, p):
        if ep == "target.json":
            return paged("targets", [], p)
        assert ep == "target/search.json" and p["q"] == "epidermal growth factor"
        return paged("targets", EGFR_TARGETS, p)

    client, session = make_client(handler)
    df = find_targets("epidermal growth factor", client=client)
    assert df.loc[0, "target_chembl_id"] == "CHEMBL203"
    assert set(df["match"]) == {"free_text"}


def test_uniprot_lookup(make_client):
    def handler(ep, p):
        assert p["target_components__accession"] == "P00533"
        return paged("targets", EGFR_TARGETS[:2] + EGFR_TARGETS[3:], p)

    client, _ = make_client(handler)
    df = find_targets("P00533", client=client)
    assert df.loc[0, "target_chembl_id"] == "CHEMBL203"
    assert set(df["match"]) == {"uniprot"}


def test_chembl_id_lookup_and_missing(make_client):
    def handler(ep, p):
        from conftest import FakeResponse
        if ep == "target/CHEMBL203.json":
            return EGFR_TARGETS[3]
        return FakeResponse(404)

    client, _ = make_client(handler)
    assert list(find_targets("chembl203", client=client)["target_chembl_id"]) == ["CHEMBL203"]
    assert find_targets("CHEMBL999999999", client=client).empty


def test_resolve_target(make_client):
    client, _ = make_client(lambda ep, p: paged("targets", EGFR_TARGETS, p))
    best = resolve_target("EGFR", client=client)
    assert best["target_chembl_id"] == "CHEMBL203"
    assert best["n_candidates"] == 4


def test_resolve_target_not_found(make_client):
    client, _ = make_client(lambda ep, p: paged("targets", [], p))
    with pytest.raises(LookupError):
        resolve_target("NOTAREALGENE", client=client)
