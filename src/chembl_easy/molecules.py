"""Find compounds by name, ChEMBL ID, similarity or substructure.

Problems handled:

* Name lookups: drug names live in ``pref_name`` *and* in synonyms/trade
  names (``Gleevec`` -> imatinib). We try exact synonym matches first, then
  ChEMBL's full-text search.
* SMILES in URLs: pasting ``C#C...`` or ``F/C=C/F`` into the URL path breaks
  the request (``#`` starts a URL fragment -> 404). Queries are always sent
  as encoded parameters.
* Invalid SMILES produce a server-side 500; we validate with RDKit first.
* Similarity thresholds below 40 % are rejected by ChEMBL; we check early.
"""

from __future__ import annotations

import pandas as pd

from chembl_easy.chem import validate_smiles
from chembl_easy.client import ChEMBLClient, NotFoundError
from chembl_easy.targets import CHEMBL_ID_RE

MOLECULE_FIELDS = (
    "molecule_chembl_id,pref_name,max_phase,molecule_type,"
    "molecule_structures,molecule_properties"
)
MOLECULE_COLUMNS = [
    "molecule_chembl_id",
    "pref_name",
    "max_phase",
    "molecule_type",
    "canonical_smiles",
    "standard_inchi_key",
    "full_mwt",
    "alogp",
    "hba",
    "hbd",
    "psa",
    "rtb",
    "num_ro5_violations",
    "qed_weighted",
]


def flatten_molecule(rec: dict) -> dict:
    """Flatten a ChEMBL molecule record into one flat row (no molfile)."""
    structs = rec.get("molecule_structures") or {}
    props = rec.get("molecule_properties") or {}
    row = {
        "molecule_chembl_id": rec.get("molecule_chembl_id"),
        "pref_name": rec.get("pref_name"),
        "max_phase": rec.get("max_phase"),
        "molecule_type": rec.get("molecule_type"),
        "canonical_smiles": structs.get("canonical_smiles"),
        "standard_inchi_key": structs.get("standard_inchi_key"),
    }
    for key in MOLECULE_COLUMNS[6:]:
        row[key] = props.get(key)
    if "similarity" in rec:
        row["similarity"] = float(rec["similarity"])
    return row


def _frame(rows: list[dict], extra: list[str] | None = None) -> pd.DataFrame:
    cols = MOLECULE_COLUMNS + (extra or [])
    df = pd.DataFrame(rows, columns=cols)
    for col in ("full_mwt", "alogp", "psa", "qed_weighted"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def get_molecule(chembl_id: str, client: ChEMBLClient | None = None) -> dict:
    """One molecule by ChEMBL ID, flattened."""
    client = client or ChEMBLClient()
    rec = client.get(f"molecule/{chembl_id.strip().upper()}", {"only": MOLECULE_FIELDS})
    return flatten_molecule(rec)


def search_molecules(query: str, client: ChEMBLClient | None = None, max_records: int = 25) -> pd.DataFrame:
    """Search by ChEMBL ID, preferred name, synonym/trade name, or free text."""
    client = client or ChEMBLClient()
    q = query.strip()
    if CHEMBL_ID_RE.match(q):
        try:
            return _frame([get_molecule(q, client)])
        except NotFoundError:
            return _frame([])
    rows: list[dict] = []
    for params in (
        {"pref_name__iexact": q},
        {"molecule_synonyms__molecule_synonym__iexact": q},
    ):
        for rec in client.paginate(
            "molecule", {**params, "only": MOLECULE_FIELDS}, key="molecules", max_records=max_records
        ):
            rows.append(flatten_molecule(rec))
    if not rows:
        for rec in client.paginate(
            "molecule/search", {"q": q, "only": MOLECULE_FIELDS}, key="molecules", max_records=max_records
        ):
            rows.append(flatten_molecule(rec))
    return _frame(rows).drop_duplicates("molecule_chembl_id").reset_index(drop=True)


def similarity_search(
    smiles: str,
    threshold: int = 70,
    client: ChEMBLClient | None = None,
    max_records: int = 100,
) -> pd.DataFrame:
    """Compounds with Tanimoto similarity >= ``threshold`` % to ``smiles``.

    ChEMBL computes similarity server-side (Morgan/ECFP-like fingerprints via
    FPSim2) and accepts thresholds between 40 and 100.
    """
    if not 40 <= int(threshold) <= 100:
        raise ValueError("ChEMBL similarity threshold must be between 40 and 100 (percent)")
    validate_smiles(smiles)
    client = client or ChEMBLClient()
    rows = [
        flatten_molecule(rec)
        for rec in client.paginate(
            "similarity",
            {"smiles": smiles.strip(), "similarity": int(threshold), "only": MOLECULE_FIELDS + ",similarity"},
            key="molecules",
            max_records=max_records,
        )
    ]
    df = _frame(rows, extra=["similarity"])
    return df.sort_values("similarity", ascending=False, kind="stable").reset_index(drop=True)


def substructure_search(
    smiles: str,
    client: ChEMBLClient | None = None,
    max_records: int = 100,
) -> pd.DataFrame:
    """Compounds containing ``smiles`` as a substructure (first ``max_records`` hits).

    Generic cores can match tens of thousands of compounds; keep
    ``max_records`` modest and refine the query instead of downloading all.
    """
    validate_smiles(smiles)
    client = client or ChEMBLClient()
    rows = [
        flatten_molecule(rec)
        for rec in client.paginate(
            "substructure",
            {"smiles": smiles.strip(), "only": MOLECULE_FIELDS},
            key="molecules",
            max_records=max_records,
        )
    ]
    return _frame(rows)
