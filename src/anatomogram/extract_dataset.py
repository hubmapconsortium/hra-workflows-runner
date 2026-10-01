import argparse
import json
import os
import re
from pathlib import Path

import anndata
import pandas as pd

from anatomogram_utils import INDIVIDUAL_COLUMN, get_dataset_ids

COUNTS_LAYER = os.environ.get("ANATOMOGRAM_COUNTS_LAYER", "filtered")
GENE_SYMBOL_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_GENE_SYMBOL", "gene_symbols")
ORGAN_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_ORGAN", "organism_part")
SEX_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_SEX", "sex")
AGE_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_AGE", "age")
RACE_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_RACE", "ethnic_group")
DISEASE_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_DISEASE", "disease")


def _get_value(obs: pd.DataFrame, column: str):
    if column not in obs.columns or obs[column].isna().all():
        return None
    return str(obs[column].dropna().iloc[0])


def _parse_age(value):
    match = re.match(r"^\s*(\d+(?:\.\d+)?)\s*year", value or "")
    return float(match[1]) if match else None


def main(args: argparse.Namespace):
    """Subsets and prints metadata from a h5ad file.

    Args:
        args (argparse.Namespace): CLI arguments, must include "file", "experiment", "dataset", and "output"
    """
    data = anndata.read_h5ad(args.file)
    mask = get_dataset_ids(data.obs, args.experiment) == args.dataset
    subset = data[mask].copy()

    # Drop precomputed embeddings, graphs, and markers
    subset.obsm.clear()
    subset.obsp.clear()
    subset.varm.clear()
    subset.uns.clear()

    # X is scaled data. Replace it with the filtered counts and drop the other layers
    counts = subset.layers[COUNTS_LAYER]
    subset.layers.clear()
    subset.X = counts
    subset.layers["counts"] = counts
    subset.var = subset.var.rename(columns={GENE_SYMBOL_COLUMN: "feature_name"})
    subset.write_h5ad(args.output, compression="gzip")

    obs = subset.obs
    age = _parse_age(_get_value(obs, AGE_COLUMN))
    metadata = {
        "organ": _get_value(obs, ORGAN_COLUMN),
        "sex": _get_value(obs, SEX_COLUMN),
        "age": int(age) if age is not None and age.is_integer() else age,
        "race": _get_value(obs, RACE_COLUMN),
        "disease": _get_value(obs, DISEASE_COLUMN),
        "donor_id": _get_value(obs, INDIVIDUAL_COLUMN),
        "cell_count": len(obs),
        "gene_count": len(subset.var),
    }
    print(json.dumps(metadata, indent=2))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract an anatomogram dataset from an experiment h5ad file")
    parser.add_argument("file", type=Path, help="Experiment h5ad file")
    parser.add_argument("--experiment", required=True, help="Experiment accession code")
    parser.add_argument("--dataset", required=True, help="Dataset id")
    parser.add_argument("--output", type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
