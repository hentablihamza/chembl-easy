"""Download and quality-filter bioactivity records for one target.

Problems handled here:

* Filtering on the **server** (target, activity type, relation, assay type,
  ``pchembl_value`` present) so you download thousands, not hundreds of
  thousands, of rows.
* Requesting only the useful columns (``only=``) - full activity records are
  large and slow.
* **Assay confidence score** lives on the *assay*, not the activity: we fetch
  the target's assays once and join, so you can keep e.g. only score 9
  (direct single-protein assignment).
* Quiet data-quality traps: mutant assays (``assay_variant_mutation``),
  records flagged by ChEMBL curators (``data_validity_comment``) and
  ``potential_duplicate`` records (the same value re-reported in a later paper).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Callable, Sequence

import pandas as pd

from chembl_easy.client import ChEMBLClient

logger = logging.getLogger(__name__)

ACTIVITY_FIELDS = [
    "activity_id",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "molecule_pref_name",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "activity_comment",
    "data_validity_comment",
    "potential_duplicate",
    "assay_chembl_id",
    "assay_type",
    "assay_description",
    "assay_variant_mutation",
    "target_chembl_id",
    "target_organism",
    "document_chembl_id",
    "document_year",
]


@dataclass
class ActivityQuery:
    """Filters for :func:`fetch_activities`. Defaults suit a typical QSAR dataset."""

    standard_types: Sequence[str] = ("IC50",)
    relations: Sequence[str] | None = ("=",)
    assay_types: Sequence[str] | None = ("B",)
    pchembl_only: bool = False
    min_confidence: int | None = None
    exclude_mutants: bool = True
    exclude_flagged: bool = True
    exclude_potential_duplicates: bool = True

    def server_params(self, target_chembl_id: str) -> dict:
        params: dict = {"target_chembl_id": target_chembl_id, "only": ",".join(ACTIVITY_FIELDS)}
        if self.standard_types:
            params["standard_type__in"] = ",".join(self.standard_types)
        if self.relations:
            params["standard_relation__in"] = ",".join(self.relations)
        if self.assay_types:
            params["assay_type__in"] = ",".join(self.assay_types)
        if self.pchembl_only:
            params["pchembl_value__isnull"] = "false"
        return params


@dataclass
class FetchReport:
    """Row counts at each filtering step, for the methods section of your paper."""

    target_chembl_id: str
    steps: list[tuple[str, int]] = field(default_factory=list)

    def add(self, name: str, n: int) -> None:
        self.steps.append((name, int(n)))
        logger.info("%-45s %d", name, n)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.steps, columns=["step", "rows"])

    def to_markdown(self) -> str:
        lines = ["| Step | Rows |", "| --- | ---: |"]
        lines += [f"| {name} | {n:,} |" for name, n in self.steps]
        return "\n".join(lines)


def fetch_assay_confidence(target_chembl_id: str, client: ChEMBLClient) -> pd.DataFrame:
    """Assay ChEMBL IDs and confidence scores for every assay on a target."""
    rows = list(
        client.paginate(
            "assay",
            {
                "target_chembl_id": target_chembl_id,
                "only": "assay_chembl_id,confidence_score,assay_organism",
            },
            key="assays",
        )
    )
    df = pd.DataFrame(rows, columns=["assay_chembl_id", "confidence_score", "assay_organism"])
    df["confidence_score"] = pd.to_numeric(df["confidence_score"], errors="coerce")
    return df


def fetch_activities(
    target_chembl_id: str,
    query: ActivityQuery | None = None,
    client: ChEMBLClient | None = None,
    max_records: int | None = None,
    progress: Callable[[int, int], None] | None = None,
    report: FetchReport | None = None,
) -> pd.DataFrame:
    """Download activities for ``target_chembl_id`` and apply quality filters.

    Returns one row per activity record with ChEMBL's original columns plus
    ``confidence_score`` (when ``min_confidence`` is set). Use
    :func:`chembl_easy.units.harmonize_units` and
    :func:`chembl_easy.aggregate.aggregate_by_compound` next.
    """
    client = client or ChEMBLClient()
    query = query or ActivityQuery()
    report = report if report is not None else FetchReport(target_chembl_id)

    params = query.server_params(target_chembl_id)
    total = client.count("activity", params)
    report.add("Server-side match (type/relation/assay type)", total)
    records = list(
        client.paginate("activity", params, key="activities", max_records=max_records, progress=progress)
    )
    df = pd.DataFrame(records)
    for col in ACTIVITY_FIELDS:
        if col not in df:
            df[col] = pd.Series(dtype=object)
    df = df[ACTIVITY_FIELDS]
    report.add("Downloaded", len(df))

    if query.exclude_flagged:
        df = df[df["data_validity_comment"].isna()]
        report.add("Without curator data-validity flags", len(df))
    if query.exclude_potential_duplicates:
        dup = pd.to_numeric(df["potential_duplicate"], errors="coerce").fillna(0)
        df = df[dup != 1]
        report.add("Without potential duplicates", len(df))
    if query.exclude_mutants:
        df = df[df["assay_variant_mutation"].isna()]
        report.add("Wild-type assays only (no variant/mutant)", len(df))
    if query.min_confidence is not None:
        conf = fetch_assay_confidence(target_chembl_id, client)
        df = df.merge(conf[["assay_chembl_id", "confidence_score"]], on="assay_chembl_id", how="left")
        df = df[df["confidence_score"] >= query.min_confidence]
        report.add(f"Assay confidence score >= {query.min_confidence}", len(df))
    df = df[df["canonical_smiles"].notna()]
    report.add("With a structure (SMILES)", len(df))
    return df.reset_index(drop=True)
