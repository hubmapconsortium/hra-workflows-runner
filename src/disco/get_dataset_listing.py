import argparse
import json
from pathlib import Path


def main(args: argparse.Namespace):
    """Prints all sample_ids (keys) from DISCO's map.json file."""
    with open(args.file, "r") as f:
        disco_map = json.load(f)
    
    # Get all keys (which are sample IDs)
    sample_ids = list(disco_map.keys())

    print(json.dumps(sample_ids))


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Get the dataset listing for DISCO")
    parser.add_argument("file", type=Path, help="Path to map.json file")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
