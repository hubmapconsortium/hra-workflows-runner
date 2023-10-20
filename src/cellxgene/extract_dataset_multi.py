import argparse
import contextlib
import json
import logging
import os
import re
import typing as t
from pathlib import Path

import anndata
from anndata import AnnData


StrList = t.List[str]

class Asset(t.TypedDict):
    id: str
    dataset: str


class Dataset(t.TypedDict):
    id: str
    donor: str
    tissue: str
    outputFile: str


class ExtractInfo(t.TypedDict):
    assets: t.List[Asset]
    assetFiles: t.List[str]
    datasets: t.List[Dataset]


_logger = logging.getLogger(__name__)


@contextlib.contextmanager
def _log_event(format: str, *args):
    try:
        _logger.info(format, "Start", *args)
        yield None
        _logger.info(format, "End", *args)
    except Exception as error:
        _logger.error(format, "Failure", *args, exc_info=error)
        raise


def _get_env_list(key: str, default: StrList) -> StrList:
    raw = os.environ.get(key)
    return re.split("[\s,;]+", raw) if raw else default


DONOR_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_DONOR", ["donor_id"])
TISSUE_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_TISSUE", ["tissue"])


class AssetSplitter:
    def __init__(self, id: str, file: Path, *column_candidates: StrList):
        self.id = id
        self.file = file
        self.column_candidates = column_candidates

    def split(self, splits: t.Iterable[StrList], out_dir: Path):
        with _log_event("CellXGene:AssetSplitter:%s - %s", self.id):
            with _log_event("CellXGene:AssetSplitter:DataLoading:%s - %s", self.file):
                matrix = anndata.read_h5ad(self.file)
            columns = self.get_columns(matrix)
            return [
                self.__split_single(matrix, columns, values, out_dir)
                for values in splits
            ]

    def __split_single(self, matrix: AnnData, columns: StrList, values: StrList, out_dir: Path):
        with _log_event(f"CellXGene:AssetSplitter:Splitting:{':'.join(values)}:%s"):
            file_name = "-".join(values) + ".h5ad"
            out_path = out_dir / self.id / file_name
            if out_path.exists():
                return out_path
            
            mask = self.create_mask(matrix, columns, values)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                matrix[mask].write_h5ad(out_path)
                return out_path
            except:
                # Don't leave partially written files
                out_path.unlink(missing_ok=True)
                raise

    def get_columns(self, matrix: AnnData):
        return [
            self.find_column(matrix, candidates)
            for candidates in self.column_candidates
        ]

    def find_column(self, matrix: AnnData, candidates: StrList):
        for column in candidates:
            if column in matrix.obs.columns:
                return column

        msg = f"Missing column - Candidates are: {candidates}"
        raise ValueError(msg)
    
    def create_mask(self, matrix: AnnData, columns: StrList, values: StrList):
        result = None
        for column, value in zip(columns, values):
            mask = matrix.obs[column] == value
            result = result & mask if result is not None else mask
        return result
    

class AssetCombiner:
    def __init__(self, id: str, out_file: Path, sources: StrList, files: t.List[Path]):
        self.id = id
        self.out_file = out_file
        self.sources = sources
        self.files = files

    def combine(self):
        with _log_event("CellXGene:AssetCombiner:%s - %s", self.id):
            matrices = [anndata.read_h5ad(file) for file in self.files]
            combined = anndata.concat(matrices, join="outer", label="source", keys=self.sources)
            self.out_file.parent.mkdir(parents=True, exist_ok=True)
            combined.write_h5ad(self.out_file)


def main(args: argparse.Namespace):
    with open(args.info) as info_file:
        info: ExtractInfo = json.load(info_file)

    assets = info["assets"]
    asset_files = info["assetFiles"]
    datasets = info["datasets"]

    splitters = [
        AssetSplitter(asset["id"], Path(file), DONOR_COLUMNS, TISSUE_COLUMNS)
        for asset, file in zip(assets, asset_files)
    ]
    splits = [
        (dataset["donor"], dataset["tissue"])
        for dataset in datasets
    ]

    with _log_event("CellXGene:Splitting:%s - %d source h5ad times %d splits", len(splitters), len(splits)):
        files = [
            splitter.split(splits, args.tmp_dir)
            for splitter in splitters
        ]
    
    sources = [asset["id"] for asset in assets]
    combiners = [
        AssetCombiner(dataset["id"], Path(dataset["outputFile"]), sources, _get_files(files, index))
        for index, dataset in enumerate(datasets)
    ]

    with _log_event("CellXGene:Combining:%s - %d output files", len(combiners)):
        for combiner in combiners:
            combiner.combine()


def _get_files(files: t.List[t.List[Path]], index: int):
    return [row[index] for row in files]


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
    parser.add_argument(
        "--log-level",
        type=int,
        default=logging.ERROR,
        help="Log level"
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)
    logging.captureWarnings(True)
    main(args)
