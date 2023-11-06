# https://github.com/hubmapconsortium/expr-h5ad-adjust/blob/main/bin/expr_h5ad_adjust.py

import argparse
import anndata
from pathlib import Path
from enum import Enum


class AnnDataLayer(str, Enum):
    SPLICED = "spliced"
    UNSPLICED = "unspliced"
    SPLICED_UNSPLICED_SUM = "spliced_unspliced_sum"


ASSAY_TO_LAYER_MAP = {
    "salmon_sn_rnaseq_10x": AnnDataLayer.SPLICED_UNSPLICED_SUM,
    "salmon_rnaseq_snareseq": AnnDataLayer.SPLICED_UNSPLICED_SUM,
    "salmon_rnaseq_slideseq": AnnDataLayer.SPLICED_UNSPLICED_SUM,
    "salmon_rnaseq_10x": AnnDataLayer.SPLICED,
    "salmon_rnaseq_sciseq": AnnDataLayer.SPLICED_UNSPLICED_SUM,
}


def main(args: argparse.Namespace):
    layer = ASSAY_TO_LAYER_MAP[args.assay]
    adata = anndata.read_h5ad(args.file)

    if layer in adata.layers:
        print("Replacing AnnData.X with layer", layer)
        adata.X = adata.layers[layer]
    else:
        raise ValueError(f"Layer {layer} not found")
    adata.write_h5ad(args.output)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="adjust the expr.h5ad X matrix")
    parser.add_argument("file", type=Path, help="Main data h5ad file")
    parser.add_argument("--assay", type=str, help="Assay type")
    parser.add_argument("--output", type=Path, help="Output file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
