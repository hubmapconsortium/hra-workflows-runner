import argparse
import anndata

TISSUE_COLUMN = "Tissue"
MAPPING = {
    "Bladder": "UBERON:0001255",
    "Blood": "UBERON:0000178",
    "Bone_Marrow": "UBERON:0002371",
    "Eye": "UBERON:0000970",
    "Heart": "UBERON:0000948",
    "Large_Intestine": "UBERON:0000059",
    "Liver": "UBERON:0002107",
    "Lung": "UBERON:0002048",
    "Lymph_Node": "UBERON:0000029",
    "Mammary": "UBERON:0001911",
    "Pancreas": "UBERON:0001264",
    "Prostate": "UBERON:0002367",
    "Skin": "UBERON:0002097",
    "Small_Intestine": "UBERON:0002108",
    "Spleen": "UBERON:0002106",
    "Thymus": "UBERON:0002370",
    "Trachea": "UBERON:0003126",
    "Uterus": "UBERON:0000995",
    "Vasculature": "UBERON:0004537",
}


def main(args: argparse.Namespace):
    organ = args.file.obs[TISSUE_COLUMN][0]
    print(f"organ: {MAPPING[organ]}")


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract gtex organ from a single h5ad file"
    )
    parser.add_argument("file", type=anndata.read_h5ad, help="Main data h5ad file")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
