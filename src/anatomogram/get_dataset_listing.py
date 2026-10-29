import argparse
import json
from pathlib import Path

import anndata

from anatomogram_utils import get_dataset_ids


def main(args: argparse.Namespace):
    """Prints unique dataset ids from a h5ad.

    Args:
        args (argparse.Namespace): CLI arguments, must contain "file" and "experiment"
    """
    data = anndata.read_h5ad(args.file, backed="r")
    datasets = get_dataset_ids(data.obs, args.experiment).drop_duplicates().tolist()
    print(json.dumps(datasets))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Get the dataset listing for an anatomogram experiment")
    parser.add_argument("file", type=Path, help="Experiment h5ad file")
    parser.add_argument("--experiment", required=True, help="Experiment accession code")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
