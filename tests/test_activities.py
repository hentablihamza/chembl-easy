import pytest

from chembl_easy.activities import ActivityQuery, FetchReport, fetch_activities

from conftest import activity, paged

RECORDS = [
    activity(1),
    activity(2, data_validity_comment="Outside typical range"),
    activity(3, potential_duplicate=1),
    activity(4, assay_variant_mutation="T790M"),
    activity(5, assay_chembl_id="CHEMBL_A8"),  # confidence 8
    activity(6, canonical_smiles=None),
    activity(7),
]
ASSAYS = [
    {"assay_chembl_id": "CHEMBL_A1", "confidence_score": 9, "assay_organism": "Homo sapiens"},
    {"assay_chembl_id": "CHEMBL_A8", "confidence_score": 8, "assay_organism": "Homo sapiens"},
]


def handler(ep, p):
    if ep == "activity.json":
        return paged("activities", RECORDS, p)
    if ep == "assay.json":
        return paged("assays", ASSAYS, p)
    raise AssertionError(ep)


def test_server_params():
    q = ActivityQuery(standard_types=["IC50", "Ki"], relations=["="], assay_types=["B", "F"], pchembl_only=True)
    p = q.server_params("CHEMBL203")
    assert p["target_chembl_id"] == "CHEMBL203"
    assert p["standard_type__in"] == "IC50,Ki"
    assert p["standard_relation__in"] == "="
    assert p["assay_type__in"] == "B,F"
    assert p["pchembl_value__isnull"] == "false"
    assert "canonical_smiles" in p["only"]


def test_server_params_optional_filters_absent():
    p = ActivityQuery(relations=None, assay_types=None).server_params("CHEMBL1")
    assert "standard_relation__in" not in p and "assay_type__in" not in p and "pchembl_value__isnull" not in p


def test_default_quality_filters(make_client):
    client, _ = make_client(handler)
    report = FetchReport("CHEMBL203")
    df = fetch_activities("CHEMBL203", client=client, report=report)
    assert sorted(df["activity_id"]) == [1, 5, 7]
    steps = dict(report.steps)
    assert steps["Downloaded"] == 7
    assert steps["Without curator data-validity flags"] == 6
    assert steps["Without potential duplicates"] == 5
    assert steps["Wild-type assays only (no variant/mutant)"] == 4


def test_confidence_filter(make_client):
    client, _ = make_client(handler)
    df = fetch_activities("CHEMBL203", query=ActivityQuery(min_confidence=9), client=client)
    assert sorted(df["activity_id"]) == [1, 7]
    assert (df["confidence_score"] == 9).all()


def test_keep_everything(make_client):
    client, _ = make_client(handler)
    q = ActivityQuery(exclude_mutants=False, exclude_flagged=False, exclude_potential_duplicates=False)
    df = fetch_activities("CHEMBL203", query=q, client=client)
    assert len(df) == 6  # only the record without SMILES is dropped


def test_max_records(make_client):
    client, _ = make_client(handler)
    df = fetch_activities("CHEMBL203", client=client, max_records=2)
    assert list(df["activity_id"]) == [1]


def test_empty_result_has_columns(make_client):
    client, _ = make_client(lambda ep, p: paged("activities", [], p))
    df = fetch_activities("CHEMBL1", client=client)
    assert df.empty and "canonical_smiles" in df.columns


def test_report_markdown():
    r = FetchReport("CHEMBL203")
    r.add("Downloaded", 18454)
    md = r.to_markdown()
    assert "| Downloaded | 18,454 |" in md
    assert list(r.to_frame().columns) == ["step", "rows"]
