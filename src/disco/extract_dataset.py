import json
import os
import pandas as pd
import argparse
import scanpy as sc

def main(args):
    append_meta_and_convert_to_h5ad(args.dataset, args.metadata, args.output)

def append_meta_and_convert_to_h5ad(extracted_h5_path, args_metadata, args_output):
    adata = sc.read_10x_h5(extracted_h5_path)
    meta_df = pd.read_csv(args_metadata, sep='\t', dtype=str)
    sample_id = os.path.basename(extracted_h5_path).replace('.h5', '')
    matched = meta_df[meta_df['sample_id'] == sample_id]
    # Append metadata to adata.obs if sample_id matches
    for col in matched.columns:
        adata.obs[col] = matched.iloc[0][col]
    # Drop null values
    adata.obs = adata.obs.dropna(axis=1, how='all')
    # Rename gene_ids to feature_name to be consistent with other h5ad datasets
    adata.var.rename(columns={'gene_ids': 'feature_name'}, inplace=True) 
    # Write to .h5ad
    adata.write(args_output)
    
    # Prepare summary
    summary = {
        "cell_count": adata.n_obs,
        "gene_count": adata.n_vars
    }
    
    print(json.dumps(summary, indent=2))
    return summary

def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Extract DISCO dataset from .tar.gz batches")
    parser.add_argument("--metadata", required=True, type=str, help="Path to metadata TSV file")
    parser.add_argument("--dataset", required=True, help="Path to extracted .h5 file")
    parser.add_argument("--output", required=True, type=str, help="Output directory")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)