"""Command-line interface: ``chembl-easy <command> ...``."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

from chembl_easy.activities import ActivityQuery
from chembl_easy.client import DEFAULT_CACHE_DIR, ChEMBLClient, ChEMBLError
from chembl_easy.chem import InvalidSmilesError
from chembl_easy.export import ml_ready, save_table
from chembl_easy.molecules import search_molecules, similarity_search, substructure_search
from chembl_easy.pipeline import build_target_dataset
from chembl_easy.targets import find_targets
from chembl_easy.version import __version__


def _client(args: argparse.Namespace) -> ChEMBLClient:
    cache = None if args.no_cache else Path(args.cache_dir)
    return ChEMBLClient(cache_dir=cache, timeout=args.timeout, max_retries=args.retries)


def _show(df: pd.DataFrame, out: str | None, columns: list[str] | None = None) -> None:
    if out:
        path = save_table(df, out)
        print(f"Wrote {len(df):,} rows to {path}")
        return
    view = df[[c for c in (columns or df.columns) if c in df.columns]]
    with pd.option_context("display.max_rows", 50, "display.width", 200, "display.max_colwidth", 60):
        print(view.to_string(index=False) if len(view) else "(no results)")


def _progress(done: int, total: int) -> None:
    print(f"\r  downloaded {done:,}/{total:,} records", end="", file=sys.stderr, flush=True)
    if done >= total:
        print(file=sys.stderr)


def cmd_target(args: argparse.Namespace) -> int:
    df = find_targets(args.query, client=_client(args), organism=args.organism,
                      single_protein_only=args.single_protein)
    _show(df, args.output)
    return 0 if len(df) else 1


def cmd_activities(args: argparse.Namespace) -> int:
    query = ActivityQuery(
        standard_types=args.types,
        relations=None if args.any_relation else ["="],
        assay_types=args.assay_types,
        pchembl_only=args.pchembl_only,
        min_confidence=args.min_confidence,
        exclude_mutants=not args.keep_mutants,
        exclude_flagged=not args.keep_flagged,
        exclude_potential_duplicates=not args.keep_duplicates,
    )
    ds = build_target_dataset(
        args.target, query=query, client=_client(args), organism=args.organism,
        max_range=args.max_range, drop_high_variance=args.drop_high_variance,
        max_records=args.max_records, progress=None if args.quiet else _progress,
    )
    t = ds.target
    print(f"Target: {t['target_chembl_id']} | {t['pref_name']} | {t['target_type']} | {t['organism']}"
          f" ({t['n_candidates']} candidate target(s) for {args.target!r})")
    print(ds.report.to_markdown())
    clean = ds.clean if not args.ml_ready else ml_ready(ds.clean, threshold=args.active_threshold)
    _show(clean, args.output, ["compound_id", "smiles", "p_activity", "n_measurements", "high_variance"])
    if args.raw_output:
        path = save_table(ds.raw, args.raw_output)
        print(f"Wrote {len(ds.raw):,} raw activity records to {path}")
    if args.summary:
        Path(args.summary).write_text(json.dumps(
            {**ds.summary(), "funnel": ds.report.steps}, indent=2), encoding="utf-8")
        print(f"Wrote summary to {args.summary}")
    return 0


def cmd_molecule(args: argparse.Namespace) -> int:
    df = search_molecules(args.query, client=_client(args), max_records=args.max_records)
    _show(df, args.output, ["molecule_chembl_id", "pref_name", "max_phase", "full_mwt", "canonical_smiles"])
    return 0 if len(df) else 1


def cmd_similar(args: argparse.Namespace) -> int:
    df = similarity_search(args.smiles, threshold=args.threshold, client=_client(args),
                           max_records=args.max_records)
    _show(df, args.output, ["molecule_chembl_id", "pref_name", "similarity", "max_phase", "canonical_smiles"])
    return 0


def cmd_substructure(args: argparse.Namespace) -> int:
    df = substructure_search(args.smiles, client=_client(args), max_records=args.max_records)
    _show(df, args.output, ["molecule_chembl_id", "pref_name", "max_phase", "canonical_smiles"])
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    print(json.dumps(_client(args).status(), indent=2))
    return 0


def cmd_clear_cache(args: argparse.Namespace) -> int:
    n = ChEMBLClient(cache_dir=Path(args.cache_dir)).clear_cache()
    print(f"Removed {n} cached responses from {args.cache_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="chembl-easy",
        description="Search ChEMBL and download clean, ML-ready bioactivity data.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR), help="on-disk response cache (default: %(default)s)")
    p.add_argument("--no-cache", action="store_true", help="always hit the network")
    p.add_argument("--timeout", type=float, default=90.0, help="per-request timeout in seconds")
    p.add_argument("--retries", type=int, default=5, help="retries for timeouts / 429 / 5xx")
    p.add_argument("-v", "--verbose", action="store_true", help="log every filtering step and retry")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("target", help="find targets by gene symbol, UniProt accession, ChEMBL ID or text")
    s.add_argument("query")
    s.add_argument("--organism", default="Homo sapiens", help="ranked first (default: %(default)s)")
    s.add_argument("--single-protein", action="store_true", help="only SINGLE PROTEIN targets")
    s.add_argument("-o", "--output")
    s.set_defaults(func=cmd_target)

    s = sub.add_parser("activities", help="download, filter, harmonize and aggregate bioactivities")
    s.add_argument("target", help="gene symbol, UniProt accession, ChEMBL target ID or name")
    s.add_argument("--types", nargs="+", default=["IC50"], help="standard types (default: IC50)")
    s.add_argument("--assay-types", nargs="+", default=["B"], help="B=binding, F=functional, A=ADMET ... (default: B)")
    s.add_argument("--any-relation", action="store_true", help="keep '<', '>' ... (default keeps '=' only)")
    s.add_argument("--pchembl-only", action="store_true", help="only records with a pChEMBL value")
    s.add_argument("--min-confidence", type=int, help="minimum assay confidence score (9 = direct single protein)")
    s.add_argument("--keep-mutants", action="store_true")
    s.add_argument("--keep-flagged", action="store_true", help="keep records with a data_validity_comment")
    s.add_argument("--keep-duplicates", action="store_true", help="keep potential_duplicate records")
    s.add_argument("--organism", default="Homo sapiens")
    s.add_argument("--max-range", type=float, default=1.0, help="flag compounds whose measurements span more log units")
    s.add_argument("--drop-high-variance", action="store_true")
    s.add_argument("--max-records", type=int, help="stop after N activity records (quick tests)")
    s.add_argument("--ml-ready", action="store_true", help="write only compound_id, smiles, p_activity[, active]")
    s.add_argument("--active-threshold", type=float, help="with --ml-ready: add active = p_activity >= threshold")
    s.add_argument("-o", "--output", help="per-compound table (.csv, .tsv or .parquet)")
    s.add_argument("--raw-output", help="also save filtered raw activity records")
    s.add_argument("--summary", help="write a JSON summary (counts, funnel, ChEMBL release)")
    s.add_argument("-q", "--quiet", action="store_true", help="no progress bar")
    s.set_defaults(func=cmd_activities)

    s = sub.add_parser("molecule", help="look up compounds by name, synonym, trade name or ChEMBL ID")
    s.add_argument("query")
    s.add_argument("--max-records", type=int, default=25)
    s.add_argument("-o", "--output")
    s.set_defaults(func=cmd_molecule)

    s = sub.add_parser("similar", help="Tanimoto similarity search (threshold 40-100 %%)")
    s.add_argument("smiles")
    s.add_argument("--threshold", type=int, default=70)
    s.add_argument("--max-records", type=int, default=100)
    s.add_argument("-o", "--output")
    s.set_defaults(func=cmd_similar)

    s = sub.add_parser("substructure", help="substructure search")
    s.add_argument("smiles")
    s.add_argument("--max-records", type=int, default=100)
    s.add_argument("-o", "--output")
    s.set_defaults(func=cmd_substructure)

    s = sub.add_parser("status", help="ChEMBL release and database counts")
    s.set_defaults(func=cmd_status)

    s = sub.add_parser("clear-cache", help="delete cached API responses")
    s.set_defaults(func=cmd_clear_cache)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    try:
        return args.func(args)
    except (InvalidSmilesError, ValueError, LookupError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ChEMBLError as exc:
        print(f"ChEMBL error: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
