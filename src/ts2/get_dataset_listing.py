import argparse
import json
import os
from pathlib import Path

import anndata

SAMPLE_ID_COLUMN = os.environ.get("TS2_COLUMN_SAMPLE_ID", "sample_id")


def main(args: argparse.Namespace):
    """Prints unique dataset ids from a h5ad.

    Args:
        args (argparse.Namespace): CLI arguments, must contain "file"
    """
    data = anndata.read_h5ad(args.file)
    datasets = data.obs[SAMPLE_ID_COLUMN].drop_duplicates().tolist()
    print(json.dumps([f"TS2-{stub}" for stub in datasets]))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Get the dataset listing for TS2")
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
