"""Shared fixtures: a fake HTTP session so no test touches the network."""

from __future__ import annotations

import json
from typing import Any, Callable

import pytest

from chembl_easy.client import ChEMBLClient

BASE = "https://www.ebi.ac.uk/chembl/api/data"


class FakeResponse:
    def __init__(self, status_code: int = 200, payload: Any = None, text: str | None = None,
                 headers: dict | None = None, url: str = ""):
        self.status_code = status_code
        self._payload = payload
        self.text = text if text is not None else (json.dumps(payload) if payload is not None else "")
        self.headers = headers or {}
        self.url = url

    def json(self):
        if self._payload is None:
            raise ValueError("no JSON")
        return self._payload


class FakeSession:
    """Routes ``get(url, params)`` to a handler and records every call."""

    def __init__(self, handler: Callable[[str, dict], Any]):
        self.handler = handler
        self.calls: list[tuple[str, dict]] = []

    def get(self, url, params=None, timeout=None):
        params = dict(params or {})
        self.calls.append((url, params))
        result = self.handler(url.replace(BASE + "/", ""), params)
        if isinstance(result, Exception):
            raise result
        if isinstance(result, FakeResponse):
            return result
        return FakeResponse(200, result, url=url)


def paged(key: str, records: list[dict], params: dict) -> dict:
    """Serve ``records`` with ChEMBL-style limit/offset pagination."""
    limit = int(params.get("limit", 20))
    offset = int(params.get("offset", 0))
    chunk = records[offset: offset + limit]
    has_next = offset + limit < len(records)
    return {
        key: chunk,
        "page_meta": {
            "limit": limit,
            "offset": offset,
            "total_count": len(records),
            "next": f"/chembl/api/data/x.json?offset={offset + limit}" if has_next else None,
            "previous": None,
        },
    }


@pytest.fixture
def make_client(tmp_path):
    """Factory: ``make_client(handler, cache=False)`` -> (client, session)."""

    def _make(handler, cache: bool = False, **kwargs):
        session = FakeSession(handler)
        client = ChEMBLClient(
            cache_dir=(tmp_path / "cache") if cache else None,
            session=session,
            sleep=lambda s: None,
            min_interval=0,
            **kwargs,
        )
        return client, session

    return _make


def activity(i: int, **overrides) -> dict:
    rec = {
        "activity_id": i,
        "molecule_chembl_id": f"CHEMBL{i}",
        "parent_molecule_chembl_id": f"CHEMBL{i}",
        "molecule_pref_name": None,
        "canonical_smiles": "c1ccccc1O",
        "standard_type": "IC50",
        "standard_relation": "=",
        "standard_value": "100.0",
        "standard_units": "nM",
        "pchembl_value": "7.00",
        "activity_comment": None,
        "data_validity_comment": None,
        "potential_duplicate": 0,
        "assay_chembl_id": "CHEMBL_A1",
        "assay_type": "B",
        "assay_description": "test",
        "assay_variant_mutation": None,
        "target_chembl_id": "CHEMBL203",
        "target_organism": "Homo sapiens",
        "document_chembl_id": "CHEMBL_D1",
        "document_year": 2010,
    }
    rec.update(overrides)
    return rec


def target(tid: str, name: str, ttype: str, organism: str = "Homo sapiens", accessions=("P00533",)) -> dict:
    return {
        "target_chembl_id": tid,
        "pref_name": name,
        "target_type": ttype,
        "organism": organism,
        "tax_id": 9606,
        "target_components": [{"accession": a} for a in accessions],
    }
