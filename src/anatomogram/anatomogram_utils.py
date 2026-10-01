import os
import re

import pandas as pd

DATASET_PREFIX = "ANATOMOGRAM"
INDIVIDUAL_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_INDIVIDUAL", "individual")
SAMPLING_SITE_COLUMN = os.environ.get("ANATOMOGRAM_COLUMN_SAMPLING_SITE", "sampling_site")


def _slug(value) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")


def get_dataset_ids(obs: pd.DataFrame, experiment: str) -> pd.Series:
    """Computes the dataset id for each cell.
    Datasets are split by donor and, when available, sampling site.

    Args:
        obs (pd.DataFrame): Cell metadata
        experiment (str): Experiment accession code, i.e. E-CURD-119

    Returns:
        pd.Series: Dataset id for each cell
    """
    ids = f"{DATASET_PREFIX}-{experiment}-" + obs[INDIVIDUAL_COLUMN].astype(str).map(_slug)
    if SAMPLING_SITE_COLUMN in obs.columns:
        ids = ids + "-" + obs[SAMPLING_SITE_COLUMN].astype(str).map(_slug)
    return ids
