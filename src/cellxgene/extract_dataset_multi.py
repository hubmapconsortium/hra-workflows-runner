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
import pandas


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


_FORCE = os.environ.get('FORCE', '').lower() not in ['', 'false']
_logger = logging.getLogger(__name__)


@contextlib.contextmanager
def _log_event(format: str, *args):
    """Logs the start, end, and any error from code run inside a with statement.

    The event is always passed before the additional arguments to the logger.

    Args:
        format (str): Formatting string passed to the logger
        *args (Any): Additional arguments passed to the logger
    """
    try:
        _logger.info(format, "Start", *args)
        yield None
        _logger.info(format, "End", *args)
    except Exception as error:
        _logger.error(format, "Failure", *args, exc_info=error)
        raise


def _get_env_list(key: str, default: StrList) -> StrList:
    """Parses an environment variable into a list.

    The parser splits the value using the separators ",", ";", and spaces.
    Multiple consecutive separators are treated as a single one.

    Args:
        key (str): Environment variable name
        default (StrList): Default value returned when the environment variable is not set

    Returns:
        StrList: The parsed list

    Examples:
        "a,b;c d" -> ["a", "b", "c", "d"]
        "a,;,b" -> ["a", "b"]
    """
    raw = os.environ.get(key)
    return re.split("[\s,;]+", raw) if raw else default


DONOR_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_DONOR", ["donor_id"])
TISSUE_COLUMNS = _get_env_list("CELLXGENE_COLUMNS_TISSUE", ["tissue"])


class AssetSplitter:
    """Splits an asset h5ad into multiple smaller file each containing
    the rows matching specific column value combinations.

    Attributes:
        id (str): Id of asset
        file (Path): Path to asset file
        column_candidates (t.List[StrList]): Possible column names to search for in the asset file
    """

    def __init__(self, id: str, file: Path, *column_candidates: StrList):
        self.id = id
        self.file = file
        self.column_candidates = column_candidates

    def split(self, splits: t.Iterable[StrList], out_dir: Path):
        """Splits the asset file.

        Args:
            splits (t.Iterable[StrList]): Combinations of values to mask each split from
            out_dir (Path): Output directory to store split files in

        Returns:
            t.List[t.Optional[Path]]: A path for each split referencing the smaller h5ad,
                None for splits that would result in empty h5ad files.
        """
        with _log_event("CellXGene:AssetSplitter:%s - %s", self.id):
            with _log_event("CellXGene:AssetSplitter:DataLoading:%s - %s", self.file):
                matrix = anndata.read_h5ad(self.file)
            columns = self.get_columns(matrix)
            return [
                self.__split_single(matrix, columns, values, out_dir)
                for values in splits
            ]

    def __split_single(self, matrix: AnnData, columns: StrList, values: StrList, out_dir: Path):
        """Perform a single split.

        Args:
            matrix (AnnData): In memory h5ad data
            columns (StrList): Columns to match split values on
            values (StrList): Split values to filter rows on
            out_dir (Path): Output directory

        Returns:
            t.Optional[Path]: The path to the resulting h5ad file or None if the subset is empty
        """
        with _log_event(f"CellXGene:AssetSplitter:Splitting:{':'.join(values)}:%s"):
            file_name = "-".join(values) + ".h5ad"
            out_path = out_dir / self.id / file_name
            if not _FORCE and out_path.exists():
                return out_path
            
            mask = self.create_mask(matrix, columns, values)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                subset = matrix[mask]
                if not subset:
                    return None
                subset.write_h5ad(out_path)
                return out_path
            except:
                # Don't leave partially written files
                out_path.unlink(missing_ok=True)
                raise

    def get_columns(self, matrix: AnnData):
        """Find columns to use during splitting from a set of candidates

        Args:
            matrix (AnnData): In memory h5ad data

        Returns:
            t.List[str]: The matching columns
        """
        return [
            self.find_column(matrix, candidates)
            for candidates in self.column_candidates
        ]

    def find_column(self, matrix: AnnData, candidates: StrList):
        """Find a column which exists within the data from a list of candidate columns.

        Args:
            matrix (AnnData): In memory h5ad data
            candidates (StrList): Candidate column names

        Raises:
            ValueError: If none of the candidates exists in the data object.

        Returns:
            str: The first matching column
        """
        for column in candidates:
            if column in matrix.obs.columns:
                return column

        msg = f"Missing column - Candidates are: {candidates}"
        raise ValueError(msg)
    
    def create_mask(self, matrix: AnnData, columns: StrList, values: StrList):
        """Create a mask which selects the rows for whose column values matches the provided values.

        Args:
            matrix (AnnData): In memory h5ad data
            columns (StrList): Columns to match on
            values (StrList): Values to match agains

        Returns:
            Any: A mask for subsetting the data matrix
        """
        result = None
        for column, value in zip(columns, values):
            mask = matrix.obs[column] == value
            result = result & mask if result is not None else mask
        return result
    

class AssetCombiner:
    """Read and combine multiple smaller h5ad into a single h5ad.

    Attributes:
        id (str): Dataset id
        out_file (Path): Path to write combined h5ad to
        sources (StrList): Ids of each asset from which the smaller h5ad was extracted
        files (t.List[t.Optional[Path]]): Path to each small h5ad or None when it should be skipped
    """

    def __init__(self, id: str, out_file: Path, sources: StrList, files: t.List[t.Optional[Path]]):
        self.id = id
        self.out_file = out_file
        self.sources = sources
        self.files = files

    def combine(self):
        """Combine split asset pieces into a single h5ad.
        """
        with _log_event("CellXGene:AssetCombiner:%s - %s [%s]", self.id, ', '.join(self.sources)):
            if not _FORCE and self.out_file.exists():
                return

            matrices = [anndata.read_h5ad(file) for file in self.files if file]
            non_empty_matrices = [matrix for matrix in matrices if matrix]
            if non_empty_matrices:
                self.fix_obsm_shapes(non_empty_matrices)

                filtered_matrices = self.filter_duplicates(non_empty_matrices)
                combined = anndata.concat(
                    filtered_matrices,
                    join="outer",
                    label="source",
                    keys=self.sources
                )

                if combined:
                    self.fix_obs_columns_dtype(combined)
                    self.out_file.parent.mkdir(parents=True, exist_ok=True)
                    combined.write_h5ad(self.out_file)
                    return

            _logger.warning(f"Subset {self.id} has zero rows")

    def fix_obsm_shapes(self, matrices: t.List[AnnData]):
        """Reshapes obsm columns to always have two dimensions to enable concatenation.

        Args:
          matrices (t.List[AnnData]): Matrices to update
        """
        for matrix in matrices:
            for key in matrix.obsm:
                array = matrix.obsm[key]
                if len(array.shape) < 2:
                    matrix.obsm[key] = array.reshape((-1, 1))

    def fix_obs_columns_dtype(self, matrix: AnnData):
        """Converts object and category columns to string to prevent errors when writing h5ad file.

        Args:
          matrix (AnnData): Matrix to update
        """
        for column in matrix.obs.columns:
            array = matrix.obs[column]
            if array.dtype in ("category", "object"):
                matrix.obs[column] = array.astype(str)

    def filter_duplicates(self, matrices: t.List[AnnData]) -> t.List[AnnData]:
        """Filters matrices to ensure unique rows.

        Args:
          matrices (t.List[AnnData]): Matrices to filter
        """
        seen_index = pandas.Index([])
        result: t.List[AnnData] = []
        for matrix in matrices:
            index = matrix.obs_names
            duplicates_mask = index.isin(seen_index)
            if duplicates_mask.any():
                matrix = matrix[~duplicates_mask]
            if matrix:
                result.append(matrix)
                seen_index = seen_index.union(matrix.obs_names)

        return result


def main(args: argparse.Namespace):
    """Splits and recombines assets into a h5ad for each dataset.
    
    The assets and datasets are read from an info json file.
    The info file must be in the format of ExtractInfo.

    Args:
        args (argparse.Namespace): CLI arguments, must contain "info" and "tmp_dir"
    """
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
        default=logging.INFO,
        help="Log level"
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level)
    logging.captureWarnings(True)
    main(args)
