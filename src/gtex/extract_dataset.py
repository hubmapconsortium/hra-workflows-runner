import argparse
import os
from pathlib import Path

import anndata

ID_COLUMN = os.environ.get("GTEX_COLUMN_ID", "Sample ID")
ORGAN_COLUMN = os.environ.get("GTEX_COLUMN_ORGAN", "Tissue")
SEX_COLUMN = os.environ.get("GTEX_COLUMN_SEX", "Sex")
AGE_COLUMN = os.environ.get("GTEX_COLUMN_AGE", "Age_bin")
DONOR_ID_COLUMN = os.environ.get("GTEX_COLUMN_DONOR", "Participant ID")
TISSUE_SITE_COLUMN = os.environ.get("GTEX_COLUMN_TISSUE_SITE", "Tissue Site Detail")


def main(args: argparse.Namespace):
    """Subsets and prints information from a h5ad file.
    Printed values include "organ", "sex", "age", "donor_id", and "tissue_site".

    Args:
        args (argparse.Namespace): CLI arguments, must include "file", "dataset", and "output"
    """
    data = anndata.read_h5ad(args.file)
    mask = data.obs[ID_COLUMN] == args.dataset
    subset = data[mask]
    subset.write_h5ad(args.output)
    print("organ:", subset.obs[ORGAN_COLUMN][0])
    print("sex:", subset.obs[SEX_COLUMN][0])
    print("age:", subset.obs[AGE_COLUMN][0])
    print("donor_id:", subset.obs[DONOR_ID_COLUMN][0])
    print("tissue_site:", subset.obs[TISSUE_SITE_COLUMN][0], flush=True)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract gtex datasets from a single h5ad file"
    )
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    parser.add_argument("--dataset", help="Dataset id")
    parser.add_argument("--output", type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
