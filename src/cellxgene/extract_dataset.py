import argparse
import os
import re
import typing as t
from pathlib import Path

import anndata


def _get_env_list(key: str, default: t.List[str]) -> t.List[str]:
    raw = os.environ.get(key)
    return re.split("[\s,;]+", raw) if raw else default


SAMPLE_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_SAMPLE", ["sample", "sample_uuid"])
DONOR_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_DONOR", ["donor_id"])
TISSUE_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_TISSUE", ["tissue"])


def _find_column(matrix: anndata.AnnData, columns: t.List[str]):
    for col in columns:
        if col in matrix.obs.columns:
            return col
    return None


def _match_column(matrix: anndata.AnnData, columns: t.List[str], value: str):
    column = _find_column(matrix, columns)
    return matrix.obs[column] == value if column else None


def create_mask(matrix: anndata.AnnData, donor: str, tissue: str, sample: str):
    if sample:
        sample_mask = _match_column(matrix, SAMPLE_COLUMNS, sample)
        if sample_mask is not None:
            return sample_mask

    donor_mask = _match_column(matrix, DONOR_COLUMNS, donor)
    tissue_mask = _match_column(matrix, TISSUE_COLUMNS, tissue)
    if donor_mask is not None and tissue_mask is not None:
        return donor_mask & tissue_mask
    return None


def main(args: argparse.Namespace):
    data = anndata.read_h5ad(args.file)
    mask = create_mask(data, args.donor, args.tissue, args.sample)
    if mask is None:
        raise ValueError("Could not filter data")

    data[mask].write_h5ad(args.output)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract cellxgene datasets per sample or per donor and tissue from a single h5ad dataset"
    )
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    parser.add_argument("--donor", help="Donor ID")
    parser.add_argument("--tissue", help="Tissue")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--output", type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
