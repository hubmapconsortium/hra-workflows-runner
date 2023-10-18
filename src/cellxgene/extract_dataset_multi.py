import argparse
import json
import os
import re
import sys
import typing as t
from pathlib import Path

import anndata


class Asset(t.TypedDict):
    id: str
    dataset: str


class Dataset(t.TypedDict):
    id: str
    donor: str
    tissue: str
    sample: t.Optional[str]
    tempFiles: t.List[t.Optional[Path]]
    outputFile: str


class ExtractInfo(t.TypedDict):
    assets: t.List[Asset]
    assetFiles: t.List[str]
    datasets: t.List[Dataset]


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
    with open(args.info) as info_file:
        info: ExtractInfo = json.load(info_file)

    datasets = info["datasets"]
    for dataset in datasets:
        dataset["tempFiles"] = []

    assets = info["assets"]
    asset_files = info["assetFiles"]
    for asset, asset_file in zip(assets, asset_files):
        print("CellXGene:Split:Start", asset["id"], flush=True)
        _split_asset(asset["id"], asset_file, datasets, args.tmp_dir)
        print('CellXGene:Split:End', asset["id"], flush=True)

    sources = [asset["dataset"] for asset in info["assets"]]
    for dataset in datasets:
        print("CellXGene:Combine:Start", dataset["id"], flush=True)
        _combine_assets(dataset, sources)
        print("CellXGene:Combine:End", dataset["id"], flush=True)


def _split_asset(asset_id: str, asset_file: str, datasets: t.List[Dataset], tmp_dir: Path):
    matrix = anndata.read_h5ad(asset_file)
    for dataset in datasets:
        donor = dataset["donor"]
        tissue = dataset["tissue"]
        sample = dataset.get("sample")
        output_file = tmp_dir / asset_id / f"{donor}-{tissue}.h5ad"
        if output_file.exists():
            continue

        mask = create_mask(matrix, donor, tissue, sample)
        if mask is not None:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            matrix[mask].write_h5ad(output_file)
            dataset["tempFiles"].append(output_file)
        else:
            dataset["tempFiles"].append(None)


def _combine_assets(dataset: Dataset, sources: t.List[str]):
    id = dataset["id"]
    temp_files = dataset["tempFiles"]
    output_file = Path(dataset["outputFile"])

    if not all(temp_files) and not any(temp_files):
        print(f"{id}: Could not extract - Unable to split sources", file=sys.stderr)
        return

    try:
        _unsafe_combine_assets(filter(None, temp_files), sources, output_file)
    except Exception as error:
        print(f"{id}: Could not extract - {error}", file=sys.stderr)


def _unsafe_combine_assets(files: t.List[Path], sources: t.List[str], output_file: Path):
    data = [anndata.read_h5ad(file) for file in files]
    combined = anndata.concat(data, join="outer", label="source", keys=sources)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output_file)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract cellxgene datasets per sample or per donor and tissue from a single h5ad dataset"
    )
    parser.add_argument("info", type=Path, help="Extract information json")
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        required=True,
        help="Directory to save temporary files"
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
