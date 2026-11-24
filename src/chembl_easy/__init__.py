"""chembl-easy: search and download clean ChEMBL data without the usual pain."""

from chembl_easy.version import __version__
from chembl_easy.activities import ActivityQuery, FetchReport, fetch_activities, fetch_assay_confidence
from chembl_easy.aggregate import add_activity_class, aggregate_by_compound
from chembl_easy.chem import InvalidSmilesError, parent_smiles, validate_smiles
from chembl_easy.client import (
    ChEMBLClient,
    ChEMBLError,
    InvalidQueryError,
    NotFoundError,
)
from chembl_easy.export import ml_ready, save_table
from chembl_easy.molecules import get_molecule, search_molecules, similarity_search, substructure_search
from chembl_easy.pipeline import TargetDataset, build_target_dataset
from chembl_easy.targets import classify_query, find_targets, resolve_target
from chembl_easy.units import harmonize_units, nM_to_p, p_to_nM, to_nM

__all__ = [
    "ActivityQuery",
    "ChEMBLClient",
    "ChEMBLError",
    "FetchReport",
    "InvalidQueryError",
    "InvalidSmilesError",
    "NotFoundError",
    "TargetDataset",
    "add_activity_class",
    "aggregate_by_compound",
    "build_target_dataset",
    "classify_query",
    "fetch_activities",
    "fetch_assay_confidence",
    "find_targets",
    "get_molecule",
    "harmonize_units",
    "ml_ready",
    "nM_to_p",
    "p_to_nM",
    "parent_smiles",
    "resolve_target",
    "save_table",
    "search_molecules",
    "similarity_search",
    "substructure_search",
    "to_nM",
    "validate_smiles",
    "__version__",
]
