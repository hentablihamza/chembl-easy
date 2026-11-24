"""Resolve a target from a gene symbol, UniProt accession, ChEMBL ID or free text.

The classic trap: searching "EGFR" returns ~20 ChEMBL targets - the human
single protein, the mouse orthologue, protein families, chimeras and
protein-protein-interaction complexes. Pulling activities for the wrong one
silently gives you a different (or mixed) dataset. :func:`find_targets`
returns all candidates ranked so that the most likely intended target - a
SINGLE PROTEIN in the requested organism - comes first, and
:func:`resolve_target` picks it while telling you what else exists.
"""

from __future__ import annotations

import logging
import re

import pandas as pd

from chembl_easy.client import ChEMBLClient, NotFoundError

logger = logging.getLogger(__name__)

CHEMBL_ID_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)
# Official UniProtKB accession pattern.
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

TARGET_FIELDS = "target_chembl_id,pref_name,target_type,organism,tax_id,target_components"
TARGET_COLUMNS = [
    "target_chembl_id",
    "pref_name",
    "target_type",
    "organism",
    "accessions",
    "n_components",
    "match",
]
TYPE_RANK = {
    "SINGLE PROTEIN": 0,
    "PROTEIN COMPLEX": 1,
    "PROTEIN FAMILY": 2,
    "SELECTIVITY GROUP": 3,
    "PROTEIN COMPLEX GROUP": 3,
    "CHIMERIC PROTEIN": 4,
    "PROTEIN-PROTEIN INTERACTION": 5,
}


def classify_query(query: str) -> str:
    """Return ``"chembl_id"``, ``"uniprot"`` or ``"text"`` for a target query."""
    q = query.strip()
    if CHEMBL_ID_RE.match(q):
        return "chembl_id"
    if UNIPROT_RE.match(q.upper()) and q == q.upper():
        return "uniprot"
    return "text"


def _flatten(target: dict, match: str) -> dict:
    comps = target.get("target_components") or []
    accessions = sorted({c.get("accession") for c in comps if c.get("accession")})
    return {
        "target_chembl_id": target.get("target_chembl_id"),
        "pref_name": target.get("pref_name"),
        "target_type": target.get("target_type"),
        "organism": target.get("organism"),
        "accessions": ";".join(accessions),
        "n_components": len(comps),
        "match": match,
    }


def rank_targets(df: pd.DataFrame, organism: str | None = None) -> pd.DataFrame:
    """Sort candidates: organism match, then single protein > complex > family > other."""
    if df.empty:
        return df
    out = df.copy()
    org_rank = 0
    if organism:
        org_rank = (out["organism"].fillna("").str.lower() != organism.lower()).astype(int)
    out["_org"] = org_rank
    out["_type"] = out["target_type"].map(TYPE_RANK).fillna(9)
    out["_match"] = out["match"].map(
        {"chembl_id": 0, "uniprot": 0, "gene_symbol": 1, "free_text": 2}
    ).fillna(3)
    out = out.sort_values(["_org", "_type", "_match", "n_components"], kind="stable")
    return out.drop(columns=["_org", "_type", "_match"]).reset_index(drop=True)


def find_targets(
    query: str,
    client: ChEMBLClient | None = None,
    organism: str | None = "Homo sapiens",
    single_protein_only: bool = False,
    max_records: int = 200,
) -> pd.DataFrame:
    """Return all ChEMBL targets matching ``query``, best candidate first.

    ``query`` may be a ChEMBL target ID (``CHEMBL203``), a UniProt accession
    (``P00533``), a gene symbol (``EGFR``) or free text (``"epidermal growth
    factor receptor"``). ``organism`` only affects ranking, not filtering.
    """
    client = client or ChEMBLClient()
    kind = classify_query(query)
    q = query.strip()
    rows: list[dict] = []
    if kind == "chembl_id":
        try:
            rec = client.get(f"target/{q.upper()}", {"only": TARGET_FIELDS})
        except NotFoundError:
            rec = None
        if rec:
            rows.append(_flatten(rec, "chembl_id"))
    elif kind == "uniprot":
        for rec in client.paginate(
            "target",
            {"target_components__accession": q, "only": TARGET_FIELDS},
            key="targets",
            max_records=max_records,
        ):
            rows.append(_flatten(rec, "uniprot"))
    if kind == "text" or (kind == "uniprot" and not rows):
        for rec in client.paginate(
            "target",
            {
                "target_components__target_component_synonyms__component_synonym__iexact": q,
                "only": TARGET_FIELDS,
            },
            key="targets",
            max_records=max_records,
        ):
            rows.append(_flatten(rec, "gene_symbol"))
        if not rows:
            for rec in client.paginate(
                "target/search", {"q": q, "only": TARGET_FIELDS}, key="targets", max_records=max_records
            ):
                rows.append(_flatten(rec, "free_text"))

    df = pd.DataFrame(rows, columns=TARGET_COLUMNS)
    df = df.drop_duplicates("target_chembl_id")
    if single_protein_only:
        df = df[df["target_type"] == "SINGLE PROTEIN"]
    return rank_targets(df, organism)


def resolve_target(
    query: str,
    client: ChEMBLClient | None = None,
    organism: str | None = "Homo sapiens",
) -> dict:
    """Pick the single most plausible target for ``query`` and return it as a dict.

    An explicit ChEMBL ID is always honoured. Otherwise the top-ranked
    candidate from :func:`find_targets` is returned and the alternatives are
    logged, so you notice when "EGFR" could have meant a protein family.
    """
    candidates = find_targets(query, client=client, organism=organism)
    if candidates.empty:
        raise LookupError(f"No ChEMBL target found for {query!r}")
    best = candidates.iloc[0].to_dict()
    if len(candidates) > 1:
        others = ", ".join(
            f"{r.target_chembl_id} ({r.target_type}, {r.organism})"
            for r in candidates.iloc[1:6].itertuples()
        )
        logger.info(
            "%r resolved to %s (%s, %s). %d other candidate(s): %s",
            query, best["target_chembl_id"], best["target_type"], best["organism"],
            len(candidates) - 1, others,
        )
    best["n_candidates"] = int(len(candidates))
    return best
