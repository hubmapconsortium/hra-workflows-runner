import argparse
import json
import os
from pathlib import Path
import re

import anndata

SAMPLE_ID_COLUMN = os.environ.get("TS2_COLUMN_SAMPLE_ID", "sample_id")
DONOR_COLUMN = os.environ.get("TS2_COLUMN_DONOR", "donor")
ORGAN_COLUMN = os.environ.get("TS2_COLUMN_ORGAN", "tissue")
AS_COLUMN = os.environ.get("TS2_COLUMN_AS", "anatomical_position")
SEX_COLUMN = os.environ.get("TS2_COLUMN_SEX", "sex")
AGE_COLUMN = os.environ.get("TS2_COLUMN_AGE", "age")
RACE_COLUMN = os.environ.get("TS2_COLUMN_AGE", "ethnicity")


def main(args: argparse.Namespace):
    """Subsets and prints information from a h5ad file.
    Printed values include "organ", "sex", "age", "donor_id", "cell_count", "gene_count" and "tissue_site".

    Args:
        args (argparse.Namespace): CLI arguments, must include "file", "dataset", and "output"
    """
    data = anndata.read_h5ad(args.file)
    sample_id = re.sub(r"^TS2-", "", args.dataset)
    mask = data.obs[SAMPLE_ID_COLUMN] == sample_id
    subset = data[mask]
    subset.write_h5ad(args.output)
    metadata = {
        "organ": subset.obs[ORGAN_COLUMN].iloc[0],
        "sex": subset.obs[SEX_COLUMN].iloc[0],
        "age": int(subset.obs[AGE_COLUMN].iloc[0]),
        "race": subset.obs[RACE_COLUMN].iloc[0],
        "donor_id": subset.obs[DONOR_COLUMN].iloc[0],
        "cell_count": len(subset.obs),
        "gene_count": len(subset.var),
        "tissue_site": subset.obs[AS_COLUMN].iloc[0],
        "rna_source": "cell",
    }
    print(json.dumps(metadata, indent=2))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract TS2 datasets from a single h5ad file"
    )
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    parser.add_argument("--dataset", help="Dataset id")
    parser.add_argument("--output", type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
