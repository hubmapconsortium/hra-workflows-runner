import argparse
import json
import os
import re
import sys
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
    with open(args.info) as infoFile:
        info = json.load(infoFile)
        print('info loaded!', flush=True)
    matrices = [anndata.read_h5ad(file) for file in info["assetFiles"]]
    sources = [asset["dataset"] for asset in info["assets"]]

    for dataset in info["datasets"]:
        result = _extract_single(matrices, dataset, sources)
        if result is not None:
            _save_matrix(result, Path(dataset["outputFile"]))


def _extract_single(matrices: t.List[anndata.AnnData], dataset: dict, sources: t.List[str]):
    try:
        return _unsafe_extract_single(matrices, dataset, sources)
    except Exception as error:
        print(f"{dataset['id']}: Could not extract - {error}", file=sys.stderr)


def _unsafe_extract_single(matrices: t.List[anndata.AnnData], dataset: dict, sources: t.List[str]):
    donor = dataset["donor"]
    tissue = dataset["tissue"]
    sample = dataset.get("sample")
    keys = []
    subsets = []

    for matrix, source in zip(matrices, sources):
        mask = create_mask(matrix, donor, tissue, sample)
        if mask is not None:
            subsets.append(matrix[mask])
            keys.append(source)
    
    if not subsets:
        raise ValueError('Failed to subset matrices')
    
    return anndata.concat(subsets, join="outer", label="source", keys=keys)


def _save_matrix(matrix: anndata.AnnData, filePath: Path):
    parentDir = filePath.parent
    parentDir.mkdir(parents=True, exist_ok=True)
    matrix.write_h5ad(filePath)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract cellxgene datasets per sample or per donor and tissue from a single h5ad dataset"
    )
    parser.add_argument("info", type=Path, help="Extract information json")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
