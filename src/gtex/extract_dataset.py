import argparse
import os
from pathlib import Path

import anndata

ID_COLUMN = os.environ.get("GTEX_COLUMN_ID", "Sample ID")
ORGAN_COLUMN = os.environ.get("GTEX_COLUMN_ORGAN", "Tissue")


def main(args: argparse.Namespace):
    data = anndata.read_h5ad(args.file)
    mask = data.obs[ID_COLUMN] == args.dataset
    subset = data[mask]
    subset.write_h5ad(args.output)
    print("organ:", subset.obs[ORGAN_COLUMN][0], flush=True)


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
