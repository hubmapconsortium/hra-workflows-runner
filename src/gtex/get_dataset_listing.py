import argparse
import os
import anndata
from pathlib import Path
import json

ID_COLUMN = os.environ.get("GTEX_COLUMN_ID", "Sample ID")


def main(args: argparse.Namespace):
    """Prints unique dataset ids from a h5ad.

    Args:
        args (argparse.Namespace): CLI arguments, must contain "file"
    """
    data = anndata.read_h5ad(args.file)
    print(json.dumps(data.obs[ID_COLUMN].drop_duplicates().tolist()))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Get the dataset listing for GTEx")
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
