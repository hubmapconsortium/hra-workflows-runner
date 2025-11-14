#!/usr/bin/env python3
"""
infer_prep_from_h5ad.py

Read a single .h5ad file and determine if the data is from 'cell' or 'nucleus'
using a priority-based approach:

1. First check: adata.obs['rna_source'] (if present, use directly)
2. Second check: adata.layers['spliced'] and ['unspliced'] (use unspliced fraction)
3. Third check: Fallback to median_frac_mito and median_frac_nuclear_lncRNA logic

"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

ENSEMBL_LOOKUP_DEFAULT = Path("/teradata/hilyas/hra-workflows-runner/hra-workflows/src/assets/ensemble-lookup.csv")



def add_ensembl_data(adata, lookup_path=ENSEMBL_LOOKUP_DEFAULT, ensembl_col="ensembl_id"):
    """
    Populate adata.var['feature_name'] from an Ensembl→symbol CSV without changing n_vars.
    - Dedupe lookup by Ensembl ID
    - Strip version suffixes on both sides
    - Use Series.map (1→1), NOT merge (1→many)
    """
    var = adata.var.copy()

    # 1) Ensure we have an Ensembl column in var
    ens_in_var = None
    for c in var.columns:
        if "ensembl" in c.lower():
            ens_in_var = c
            break
    if ens_in_var is None:
        var[ensembl_col] = var.index.astype(str)
    elif ens_in_var != ensembl_col:
        var = var.rename(columns={ens_in_var: ensembl_col})

    # Normalize Ensembl IDs: drop version suffix
    var[ensembl_col] = var[ensembl_col].astype(str).str.split(".").str[0]

    # 2) Load and normalize lookup
    df = pd.read_csv(lookup_path)
    lc = {c.lower(): c for c in df.columns}
    ens_lu = lc.get("ensembl_id") or next((c for c in df.columns if "ensembl" in c.lower()), None)
    sym_lu = lc.get("gene_name") or lc.get("symbol") or next(
        (c for c in df.columns if "symbol" in c.lower() or c.lower() in ("name","gene")), None
    )
    if ens_lu is None or sym_lu is None:
        raise ValueError("Lookup CSV must contain Ensembl and gene symbol/name columns.")

    df = df[[ens_lu, sym_lu]].dropna()
    df[ens_lu] = df[ens_lu].astype(str).str.split(".").str[0]

    # Keep a single symbol per Ensembl ID (first non-empty)
    df = df[df[sym_lu].astype(str).ne("")].drop_duplicates(subset=[ens_lu], keep="first")

    # 3) Build mapping and fill feature_name via map (preserves length)
    map_dict = pd.Series(df[sym_lu].values, index=df[ens_lu].values).to_dict()

    if "feature_name" not in var.columns:
        var["feature_name"] = pd.Series(index=var.index, dtype="object")

    var["feature_name"] = var["feature_name"].where(var["feature_name"].notna(),
                                                   var[ensembl_col].map(map_dict))

    # 4) Fallbacks: try other symbol columns, then index
    for alt in ("gene_name", "gene_symbol", "symbol", "SYMBOL"):
        if alt in var.columns:
            var["feature_name"] = var["feature_name"].fillna(var[alt].astype(str).replace("nan",""))

    var["feature_name"] = var["feature_name"].fillna(pd.Series(var.index.astype(str), index=var.index))

    # 5) Safety: enforce unchanged shape
    if len(var) != adata.n_vars:
        raise ValueError(f"Internal error: var length changed ({len(var)} != n_vars {adata.n_vars}). "
                         "Check lookup duplicates.")

    adata.var = var
    return adata


def _to_str(x):
    """Helper to convert bytes/None to string."""
    if x is None:
        return ""
    if isinstance(x, bytes):
        try:
            return x.decode("utf-8", "ignore")
        except Exception:
            return ""
    return str(x)


def get_counts(adata):
    """Get count matrix from layers['counts'] or X."""
    return adata.layers['counts'] if 'counts' in adata.layers else adata.X


def get_symbols_array(adata, ensemble_lookup=ENSEMBL_LOOKUP_DEFAULT):
    """
    Return a numpy array of gene symbols.
    If 'feature_name' is missing or empty, try to build it from Ensembl lookup.
    """
    need_feature_name = (
        ("feature_name" not in adata.var.columns) or
        (adata.var["feature_name"].isna().all()) or
        (adata.var["feature_name"].astype(str).str.len().eq(0).all())
    )

    if need_feature_name and ensemble_lookup:
        add_ensembl_data(adata, lookup_path=ensemble_lookup)
        # print(adata.var["feature_name"].unique(), file=sys.stderr)

    # final fallback: use var_names if still no feature_name
    if "feature_name" in adata.var:
        ser = adata.var["feature_name"]
    else:
        ser = pd.Series(adata.var_names, index=adata.var_names)

    symbols = ser.astype("object").map(_to_str).fillna("").values
    return symbols.astype("U")


def frac_from_mask(adata, mask, X):
    """Calculate fraction of counts from genes matching mask."""
    if not np.any(mask):
        return np.zeros(adata.n_obs, dtype=float)
    if sp.issparse(X):
        tot = np.asarray(X.sum(1)).ravel()
        sub = np.asarray(X[:, mask].sum(1)).ravel()
    else:
        tot = X.sum(1)
        sub = X[:, mask].sum(1)
    tot = np.maximum(tot, 1)
    out = sub / tot
    out[np.isnan(out)] = 0.0
    return out


def get_arr_sum(matrix):
    """Sum rows of matrix and return 1D numpy array (handles sparse)."""
    if matrix is None:
        return None
    if sp.issparse(matrix):
        s = np.asarray(matrix.sum(axis=1)).ravel()
    else:
        s = np.asarray(matrix.sum(axis=1)).ravel()
    return s


NUC_LNCRNAS = {
    "NEAT1","MALAT1","XIST",
}

def make_mito_mask(adata, sym_upper):
    var = adata.var

    # 1) precomputed boolean-ish columns
    for col in var.columns:
        cl = col.lower()
        if cl in ("mt", "is_mt", "mitochondrial"):
            s = var[col]
            if s.dtype == bool:
                return s.values
            s = s.astype(str).str.lower()
            return s.isin(["true","1","t","y","mt","mitochondrial"]).values

    # 2) feature types / biotypes
    for col in ("feature_type","feature_types","gene_biotype","biotype","type"):
        if col in var.columns:
            s = var[col].astype(str).str.lower()
            # match 'mt' or 'mitochond' substrings
            m = s.str.contains(r"\bmt\b|mitochond", na=False)
            if m.any():
                return m.values

    # 3) genome / chromosome
    for col in ("genome","chrom","chromosome","seqname"):
        if col in var.columns:
            s = var[col].astype(str).str.lower()
            m = s.isin(["mt","m","chrm","chrmt"])
            if m.any():
                return m.values

    # 4) fallback: symbols with 'MT-' prefix (human)
    return np.char.startswith(sym_upper, "MT-")


def make_nuclear_lnc_mask(sym_upper):
    return np.isin(sym_upper, [g.upper() for g in NUC_LNCRNAS])

# ---------------------------
# Priority 1: Check rna_source in obs
# ---------------------------
def check_rna_source(adata):
    """
    Priority 1: Check if rna_source column exists in adata.obs.
    Returns: (result, evidence) with result in {'cell','nucleus',None}
    """
    if 'rna_source' not in adata.obs.columns:
        return None, "No 'rna_source' column found in adata.obs"

    uniques = pd.unique(adata.obs['rna_source'].astype('object'))
    uniques = [str(u).lower().strip() for u in uniques if pd.notna(u)]

    if not uniques:
        return None, "rna_source column exists but contains no valid values"

    for val in uniques:
        if 'nucleus' in val or 'nuclei' in val:
            return 'nucleus', f"Found rna_source='{val}' in adata.obs - returning nucleus"

    for val in uniques:
        if 'cell' in val:
            return 'cell', f"Found rna_source='{val}' in adata.obs - returning cell"

    return None, f"rna_source values: {uniques}, but no clear 'cell' or 'nucleus' match"


# ---------------------------
# Priority 2: Check spliced/unspliced layers
# ---------------------------
def check_spliced_unspliced(adata, threshold=0.30):
    """
    High unspliced fraction (>= threshold) → nucleus, else cell.
    Returns: (result, evidence, median_frac_unspliced)
    """
    if 'unspliced' not in adata.layers or 'spliced' not in adata.layers:
        return None, "No 'spliced' and/or 'unspliced' layers found in adata.layers", None

    try:
        unsp = get_arr_sum(adata.layers['unspliced']).astype(float)
        splic = get_arr_sum(adata.layers['spliced']).astype(float)
        total = splic + unsp + 1e-9
        frac_unspliced = unsp / total
        median_frac = float(np.nanmedian(frac_unspliced))

        if median_frac >= threshold:
            return 'nucleus', f"Median unspliced fraction={median_frac:.4f} >= {threshold} → nucleus", median_frac
        else:
            return 'cell', f"Median unspliced fraction={median_frac:.4f} < {threshold} → cell", median_frac

    except Exception as e:
        return None, f"Failed to compute unspliced fraction from layers: {e}", None


# ---------------------------
# Priority 3: Fallback to median_frac_mito and median_frac_nuclear_lncRNA
# ---------------------------
def check_frac_metrics(adata):
    """
    Fallback rule (human only):
      - If median_frac_mito < 0.020 AND median_frac_nuclear_lncRNA > 0.033 → nucleus
      - Else → cell
    Returns: (result, evidence, metrics_dict) or (None, evidence, None) if no markers found.
    """
    try:
        X = get_counts(adata)
        symbols = get_symbols_array(adata)  
        sym_upper = np.char.upper(symbols)

        mito_mask = make_mito_mask(adata, sym_upper)
        nuc_mask  = make_nuclear_lnc_mask(sym_upper)

        # If we truly can't find either marker set, don't force a call
        if not np.any(mito_mask) and not np.any(nuc_mask):
            return (None,
                    "No mitochondrial or nuclear-retained lncRNA markers detected after all lookups; "
                    "skipping fraction heuristic.",
                    None)

        # Fractions (zeros if a particular mask is empty — that's fine)
        frac_mito = frac_from_mask(adata, mito_mask, X)
        frac_nuc  = frac_from_mask(adata, nuc_mask,  X)

        median_frac_mito = float(np.median(frac_mito))
        median_frac_nuclear_lncRNA = float(np.median(frac_nuc))

        metrics = {
            "median_frac_mito": median_frac_mito,
            "median_frac_nuclear_lncRNA": median_frac_nuclear_lncRNA,
            "markers_present": {
                "mito": bool(np.any(mito_mask)),
                "nuclear_lncRNA": bool(np.any(nuc_mask))
            }
        }

        if median_frac_mito < 0.020 and median_frac_nuclear_lncRNA > 0.033:
            evidence = (f"median_frac_mito={median_frac_mito:.6f} < 0.020 AND "
                        f"median_frac_nuclear_lncRNA={median_frac_nuclear_lncRNA:.6f} > 0.033 → nucleus")
            return "nucleus", evidence, metrics
        else:
            evidence = (f"Does not meet nucleus criteria "
                        f"(mito={median_frac_mito:.6f}, nuclear_lncRNA={median_frac_nuclear_lncRNA:.6f}) → cell")
            return "cell", evidence, metrics

    except Exception as e:
        return None, f"Failed to compute fraction metrics: {e}", None


# ---------------------------
# Main inference function
# ---------------------------
def analyze_h5ad(path, unspliced_threshold=0.30):
    """
    Priority order:
      1. obs['rna_source']
      2. layers['spliced'/'unspliced']
      3. mito/lncRNA fractions
    Returns a report dict.
    """
    adata = ad.read_h5ad(path, backed=False)
    evidence = []
    result = None
    method_used = None
    metrics = {}

    # 1) rna_source
    result, evidence_msg = check_rna_source(adata)
    evidence.append(f"Priority 1 (rna_source): {evidence_msg}")
    if result:
        method_used = "rna_source_column"
        verdict = result
    else:
        # 2) spliced/unspliced
        result, evidence_msg, median_frac_unspliced = check_spliced_unspliced(adata, threshold=unspliced_threshold)
        evidence.append(f"Priority 2 (spliced/unspliced): {evidence_msg}")
        if result:
            method_used = "spliced_unspliced_layers"
            verdict = result
            if median_frac_unspliced is not None:
                metrics['median_frac_unspliced'] = median_frac_unspliced
        else:
            # 3) mito/lncRNA fallback
            result, evidence_msg, frac_metrics = check_frac_metrics(adata)
            evidence.append(f"Priority 3 (fraction_metrics): {evidence_msg}")
            if result:
                method_used = "fraction_metrics_fallback"
                verdict = result
                if frac_metrics:
                    metrics.update(frac_metrics)
            else:
                verdict = "inconclusive"
                method_used = "failed_all_methods"
                evidence.append("All three methods failed to determine cell type")

    report = {
        "input_file": path,
        "verdict": verdict,
        "method_used": method_used,
        "evidence": evidence,
        "metrics": metrics,
    }
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Infer whether an AnnData (.h5ad) file is from cells or nuclei (human-only) using priority-based inference."
    )
    parser.add_argument(
        "h5ad",
        help="Path to input .h5ad file OR dataset folder name (e.g., DISCO-1823_BA24_10x). "
             "If dataset folder name, use --base-dir to specify base directory."
    )
    parser.add_argument("--output", "-o", help="Output JSON file (default: stdout)", default=None)
    parser.add_argument(
        "--base-dir"
    )
    parser.add_argument("--unspliced-threshold", type=float, default=0.30,
                       help="Threshold for unspliced fraction to call nucleus (default: 0.30)")
    parser.add_argument("--debug", action="store_true", help="Print debug info to stderr")

    args = parser.parse_args()

    # Determine if input is a dataset folder name or full path
    input_path = Path(args.h5ad)

    # Check if it's already a full path to an .h5ad file that exists
    if input_path.is_absolute() and input_path.exists() and input_path.suffix == ".h5ad":
        h5ad_path = input_path
        if args.debug:
            print(f"Using provided full path: {h5ad_path}", file=sys.stderr)
    elif input_path.is_absolute() and input_path.exists():
        potential_h5ad = input_path / "data.h5ad"
        h5ad_path = potential_h5ad if potential_h5ad.exists() else input_path
        if args.debug:
            print(f"Resolved path: {h5ad_path}", file=sys.stderr)
    else:
        # Assume it's a dataset folder name - construct full path
        base_dir = Path(args.base_dir)
        h5ad_path = base_dir / args.h5ad / "data.h5ad"
        if args.debug:
            print(f"Treated as dataset folder name: {args.h5ad}", file=sys.stderr)
            print(f"Constructed full path: {h5ad_path}", file=sys.stderr)

    # Check if file exists
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
        report = analyze_h5ad(str(h5ad_path), unspliced_threshold=args.unspliced_threshold)
    except Exception as e:
        err = {"input_file": str(h5ad_path), "error": str(e), "verdict": "error"}
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
