#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import anndata as ad
import scipy.sparse as sp


def get_arr_sum(matrix):
    """Sum rows of matrix and return a 1D numpy array (handles sparse matrices)."""
    if matrix is None:
        return None
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix.sum(axis=1)).ravel()


def analyze_h5ad(path, unspliced_threshold=0.30):
    """
    Returns a report dict with verdict, evidence, and metrics based solely on
    the median unspliced fraction (>= threshold → nucleus).
    """
    adata = ad.read_h5ad(path, backed=False)
    evidence = []
    metrics = {}

    if "unspliced" not in adata.layers or "spliced" not in adata.layers:
        return {
            "input_file": path,
            "verdict": "inconclusive",
            "evidence": [
                "No 'spliced' and/or 'unspliced' layers found in adata.layers"
            ],
            "metrics": {},
        }

    try:
        unsp = get_arr_sum(adata.layers["unspliced"]).astype(float)
        splic = get_arr_sum(adata.layers["spliced"]).astype(float)
        total = splic + unsp + 1e-9
        frac_unspliced = unsp / total
        median_frac = float(np.nanmedian(frac_unspliced))
        metrics["median_frac_unspliced"] = median_frac

        if median_frac >= unspliced_threshold:
            verdict = "nucleus"
            evidence.append(
                f"Median unspliced fraction={median_frac:.4f} >= {unspliced_threshold} → nucleus"
            )
        else:
            verdict = "cell"
            evidence.append(
                f"Median unspliced fraction={median_frac:.4f} < {unspliced_threshold} → cell"
            )
    except Exception as exc:
        return {
            "input_file": path,
            "verdict": "error",
            "evidence": [f"Failed to compute unspliced fraction: {exc}"],
            "metrics": {},
        }

    return {
        "input_file": path,
        "verdict": verdict,
        "evidence": evidence,
        "metrics": metrics,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Infer whether an AnnData (.h5ad) file is from cells or nuclei using only spliced/unspliced layers."
    )
    parser.add_argument(
        "h5ad",
        help="Path to input .h5ad file OR dataset folder name (e.g., DISCO-1823_BA24_10x). "
        "If dataset folder name, use --base-dir to specify base directory.",
    )
    parser.add_argument(
        "--output", "-o", help="Output JSON file (default: stdout)", default=None
    )
    parser.add_argument(
        "--base-dir", help="Base directory containing dataset folders", default=None
    )
    parser.add_argument(
        "--unspliced-threshold",
        type=float,
        default=0.30,
        help="Threshold for unspliced fraction to call nucleus (default: 0.30)",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print debug info to stderr"
    )

    args = parser.parse_args()

    input_path = Path(args.h5ad)
    if input_path.is_absolute() and input_path.exists():
        if input_path.is_file() and input_path.suffix == ".h5ad":
            h5ad_path = input_path
        else:
            potential = input_path / "data.h5ad"
            h5ad_path = potential if potential.exists() else input_path
        if args.debug:
            print(f"Resolved path: {h5ad_path}", file=sys.stderr)
    else:
        if not args.base_dir:
            parser.error("--base-dir is required when passing a dataset folder name")
        base_dir = Path(args.base_dir)
        h5ad_path = base_dir / args.h5ad / "data.h5ad"
        if args.debug:
            print(
                f"Treated as dataset folder name. Constructed path: {h5ad_path}",
                file=sys.stderr,
            )

    if not h5ad_path.exists():
        err_msg = f"Error: File not found: {h5ad_path}"
        if not input_path.is_absolute():
            err_msg += f"\n  (Tried to find dataset folder '{args.h5ad}' in base directory: {args.base_dir})"
        err = {"input_file": str(h5ad_path), "error": err_msg, "verdict": "error"}
        if args.output:
            with open(args.output, "w") as fh:
                json.dump(err, fh, indent=2)
        else:
            json.dump(err, sys.stdout, indent=2)
            sys.stdout.write("\n")
        sys.exit(1)

    try:
        report = analyze_h5ad(
            str(h5ad_path), unspliced_threshold=args.unspliced_threshold
        )
    except Exception as exc:
        err = {"input_file": str(h5ad_path), "error": str(exc), "verdict": "error"}
        if args.output:
            with open(args.output, "w") as fh:
                json.dump(err, fh, indent=2)
        else:
            json.dump(err, sys.stdout, indent=2)
            sys.stdout.write("\n")
        if args.debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    if args.output:
        with open(args.output, "w") as fh:
            json.dump(report, fh, indent=2)
    else:
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")


if __name__ == "__main__":
    main()
