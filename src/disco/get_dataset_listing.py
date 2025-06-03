import argparse
import pandas as pd
import json
from pathlib import Path
import math

def main(args: argparse.Namespace):
    """Main function to process DISCO data and output sample listing"""
    # Read metadata
    meta_df = pd.read_csv(args.metadata, sep="\t", dtype=str)
    
    # Process metadata
    meta_dict = {
        row['sample_id']: {k: v for k, v in {**row.to_dict(), 'sample_id': row['sample_id']}.items() 
                          if not (isinstance(v, float) and math.isnan(v))}
        for _, row in meta_df.iterrows()
    }
    print(json.dumps(sorted(meta_dict.keys())))


def _get_arg_parser():
    parser = argparse.ArgumentParser(description='Get DISCO dataset listing')
    parser.add_argument('metadata', type=Path, help='Path to metadata TSV file')
    return parser

if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
