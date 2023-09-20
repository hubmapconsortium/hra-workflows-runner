import argparse
from pathlib import Path
import anndata

SAMPLE_ID_COLUMN = 'sample'
DONOR_ID_COLUMN = 'donor_id'
TISSUE_COLUMN = 'tissue'

def main(args: argparse.Namespace):
    data = anndata.read_h5ad(args.file)
    if args.sample:
        mask = data.obs[SAMPLE_ID_COLUMN] == args.sample
    if args.donor and args.tissue:
        mask = (data.obs[TISSUE_COLUMN] == args.tissue) & (data.obs[DONOR_ID_COLUMN] == args.donor) 
    data[mask].write_h5ad(args.output)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract cellxgene datasets per sample or per donor and tissue from a single h5ad dataset"
    )
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    parser.add_argument("--donor", help="Donor ID")
    parser.add_argument("--tissue", help="Tissue")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument('--output', type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
