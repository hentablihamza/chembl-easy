"""Live checks against the real ChEMBL API.

Skipped by default (and in CI). Run with:  CHEMBL_EASY_NETWORK_TESTS=1 pytest -m network
"""

import os

import pytest

from chembl_easy import ChEMBLClient, find_targets, search_molecules, similarity_search

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(os.environ.get("CHEMBL_EASY_NETWORK_TESTS") != "1",
                       reason="set CHEMBL_EASY_NETWORK_TESTS=1 to hit the live ChEMBL API"),
]


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    return ChEMBLClient(cache_dir=tmp_path_factory.mktemp("cache"))


def test_live_egfr_resolves_to_chembl203(client):
    df = find_targets("EGFR", client=client)
    assert df.iloc[0]["target_chembl_id"] == "CHEMBL203"


def test_live_trade_name(client):
    assert "CHEMBL941" in set(search_molecules("Gleevec", client=client)["molecule_chembl_id"])


def test_live_similarity_with_triple_bond(client):
    df = similarity_search("C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", threshold=90, client=client)
    assert "CHEMBL553" in set(df["molecule_chembl_id"])  # erlotinib
