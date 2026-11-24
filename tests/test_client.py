import pytest
import requests

from chembl_easy.client import ChEMBLClient, ChEMBLError, InvalidQueryError, NotFoundError, _cache_key

from conftest import FakeResponse, paged


def test_url_building():
    c = ChEMBLClient(cache_dir=None)
    assert c.url("activity") == "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    assert c.url("/molecule/CHEMBL25/") == "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json"
    assert c.url("status.json").endswith("/status.json")


def test_get_returns_json(make_client):
    client, session = make_client(lambda ep, p: {"ok": True, "ep": ep})
    assert client.get("status") == {"ok": True, "ep": "status.json"}
    assert len(session.calls) == 1


def test_retries_on_503_then_succeeds(make_client):
    state = {"n": 0}

    def handler(ep, p):
        state["n"] += 1
        return FakeResponse(503) if state["n"] < 3 else {"ok": 1}

    sleeps = []
    client, session = make_client(handler)
    client._sleep = sleeps.append
    assert client.get("activity") == {"ok": 1}
    assert len(session.calls) == 3
    assert sleeps == [1.5, 3.0]  # exponential backoff


def test_retries_on_timeout(make_client):
    state = {"n": 0}

    def handler(ep, p):
        state["n"] += 1
        return requests.Timeout("slow") if state["n"] == 1 else {"ok": 1}

    client, session = make_client(handler)
    assert client.get("activity") == {"ok": 1}
    assert len(session.calls) == 2


def test_retry_after_header_is_honoured(make_client):
    state = {"n": 0}

    def handler(ep, p):
        state["n"] += 1
        return FakeResponse(429, headers={"Retry-After": "7"}) if state["n"] == 1 else {"ok": 1}

    sleeps = []
    client, _ = make_client(handler)
    client._sleep = sleeps.append
    client.get("activity")
    assert sleeps == [7.0]


def test_gives_up_after_max_retries(make_client):
    client, session = make_client(lambda ep, p: FakeResponse(502), max_retries=2)
    with pytest.raises(ChEMBLError, match="after 3 attempts"):
        client.get("activity")
    assert len(session.calls) == 3


def test_404_raises_not_found_without_retry(make_client):
    client, session = make_client(lambda ep, p: FakeResponse(404, text="nope"))
    with pytest.raises(NotFoundError):
        client.get("molecule/CHEMBL0")
    assert len(session.calls) == 1


def test_bad_smiles_500_is_not_retried(make_client):
    body = 'Exception: Problem standardizing the query molecule'
    client, session = make_client(lambda ep, p: FakeResponse(500, text=body))
    with pytest.raises(InvalidQueryError):
        client.get("similarity", {"smiles": "xx", "similarity": 70})
    assert len(session.calls) == 1


def test_400_raises_immediately(make_client):
    client, session = make_client(lambda ep, p: FakeResponse(400, text="Invalid Similarity Score"))
    with pytest.raises(ChEMBLError, match="HTTP 400"):
        client.get("similarity")
    assert len(session.calls) == 1


def test_disk_cache_avoids_second_request(make_client):
    client, session = make_client(lambda ep, p: {"v": 1}, cache=True)
    assert client.get("activity", {"a": 1}) == {"v": 1}
    assert client.get("activity", {"a": 1}) == {"v": 1}
    assert len(session.calls) == 1
    assert client.n_cache_hits == 1
    client.get("activity", {"a": 2})
    assert len(session.calls) == 2
    assert client.clear_cache() == 2


def test_cache_key_is_order_independent():
    assert _cache_key("u", {"a": 1, "b": 2}) == _cache_key("u", {"b": 2, "a": 1})
    assert _cache_key("u", {"a": 1}) != _cache_key("u", {"a": 2})


def test_status_is_never_cached(make_client):
    client, session = make_client(lambda ep, p: {"chembl_db_version": "ChEMBL_37"}, cache=True)
    client.status()
    client.status()
    assert len(session.calls) == 2


def test_paginate_walks_all_pages(make_client):
    records = [{"i": i} for i in range(25)]
    client, session = make_client(lambda ep, p: paged("items", records, p))
    seen = []
    out = list(client.paginate("thing", page_size=10, progress=lambda d, t: seen.append((d, t))))
    assert [r["i"] for r in out] == list(range(25))
    assert len(session.calls) == 3
    assert seen[-1] == (25, 25)
    assert [c[1]["offset"] for c in session.calls] == [0, 10, 20]


def test_paginate_max_records_and_page_cap(make_client):
    records = [{"i": i} for i in range(5000)]
    client, session = make_client(lambda ep, p: paged("items", records, p))
    out = list(client.paginate("thing", key="items", page_size=5000, max_records=1500))
    assert len(out) == 1500
    assert all(c[1]["limit"] <= 1000 for c in session.calls)


def test_paginate_empty(make_client):
    client, _ = make_client(lambda ep, p: paged("items", [], p))
    assert list(client.paginate("thing", key="items")) == []


def test_count(make_client):
    client, session = make_client(lambda ep, p: paged("activities", [{}] * 42, p))
    assert client.count("activity", {"target_chembl_id": "CHEMBL203"}) == 42
    assert session.calls[0][1]["limit"] == 1
