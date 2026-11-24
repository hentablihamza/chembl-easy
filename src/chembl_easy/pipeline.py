"""One call from a target name to a clean, ML-ready dataset."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import pandas as pd

from chembl_easy.activities import ActivityQuery, FetchReport, fetch_activities
from chembl_easy.aggregate import aggregate_by_compound
from chembl_easy.client import ChEMBLClient
from chembl_easy.targets import resolve_target
from chembl_easy.units import harmonize_units


@dataclass
class TargetDataset:
    target: dict
    raw: pd.DataFrame
    clean: pd.DataFrame
    report: FetchReport
    chembl_release: str | None = None

    def summary(self) -> dict:
        c = self.clean
        return {
            "target_chembl_id": self.target["target_chembl_id"],
            "pref_name": self.target["pref_name"],
            "target_type": self.target["target_type"],
            "organism": self.target["organism"],
            "chembl_release": self.chembl_release,
            "activity_records": int(len(self.raw)),
            "compounds": int(len(c)),
            "high_variance_compounds": int(c["high_variance"].sum()) if len(c) else 0,
            "p_activity_median": float(c["p_activity"].median()) if len(c) else None,
            "p_activity_min": float(c["p_activity"].min()) if len(c) else None,
            "p_activity_max": float(c["p_activity"].max()) if len(c) else None,
            "actives_p_ge_6": int((c["p_activity"] >= 6).sum()) if len(c) else 0,
            "actives_p_ge_7": int((c["p_activity"] >= 7).sum()) if len(c) else 0,
        }


def build_target_dataset(
    target_query: str,
    query: ActivityQuery | None = None,
    client: ChEMBLClient | None = None,
    organism: str | None = "Homo sapiens",
    max_range: float = 1.0,
    drop_high_variance: bool = False,
    max_records: int | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> TargetDataset:
    """Resolve target -> download activities -> harmonize units -> aggregate per compound."""
    client = client or ChEMBLClient()
    target = resolve_target(target_query, client=client, organism=organism)
    tid = target["target_chembl_id"]
    report = FetchReport(tid)
    raw = fetch_activities(tid, query=query, client=client, max_records=max_records,
                           progress=progress, report=report)
    raw = harmonize_units(raw)
    usable = raw[raw["p_activity"].notna()]
    report.add("Convertible to p-scale (pIC50/pKi/...)", len(usable))
    clean = aggregate_by_compound(usable, max_range=max_range, drop_high_variance=drop_high_variance)
    report.add("Unique parent compounds", len(clean))
    try:
        release = client.status().get("chembl_db_version")
    except Exception:  # status is informative only
        release = None
    return TargetDataset(target=target, raw=raw, clean=clean, report=report, chembl_release=release)
