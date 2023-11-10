import argparse
import os
from pathlib import Path
import anndata

SEX_COLUMN = os.environ.get("GTEX_COLUMN_SEX", "sex")
AGE_COLUMN = os.environ.get("GTEX_COLUMN_AGE", "development_stage")
ETHNICITY_COLUMN = os.environ.get("GTEX_COLUMN_AGE", "self_reported_ethnicity")

def main(args: argparse.Namespace):
    data = anndata.read_h5ad(args.file)
    print("sex:", data.obs[SEX_COLUMN][0])
    print("age:", data.obs[AGE_COLUMN][0])
    print("ethnicity:", data.obs[ETHNICITY_COLUMN][0])


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract donor metadata from a h5ad file"
    )
    parser.add_argument("file", type=Path, help="data h5ad file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
