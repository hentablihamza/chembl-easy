"""RDKit helpers: SMILES validation and salt stripping."""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

rdBase.DisableLog("rdApp.*")


class InvalidSmilesError(ValueError):
    """Raised when RDKit cannot parse a SMILES string."""


def validate_smiles(smiles: str) -> str:
    """Return RDKit canonical SMILES or raise :class:`InvalidSmilesError`.

    Validate *before* calling ChEMBL similarity/substructure search: the API
    answers an unparsable SMILES with an HTTP 500 traceback, not a helpful
    message.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise InvalidSmilesError("Empty SMILES")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise InvalidSmilesError(f"RDKit cannot parse SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol)


@lru_cache(maxsize=200_000)
def parent_smiles(smiles: str) -> str | None:
    """Largest organic fragment, neutralized, as canonical SMILES (salt/solvent stripped).

    Returns ``None`` for unparsable input.
    """
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception:  # pragma: no cover - RDKit edge cases
        return None
    return Chem.MolToSmiles(mol)


@lru_cache(maxsize=200_000)
def inchikey(smiles: str) -> str | None:
    """Standard InChIKey of a SMILES, or ``None`` if it cannot be computed."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    key = Chem.MolToInchiKey(mol)
    return key or None
