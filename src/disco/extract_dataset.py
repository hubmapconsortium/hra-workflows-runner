import json
import os
import pandas as pd
import argparse
from pathlib import Path
import scanpy as sc
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed

def main(args):

    dataset_id = args.dataset
    dataset_id = dataset_id.rstrip('/').split('/')[-1] 
    extracted_h5_path = parallel_search(args.batches, dataset_id, args.output)
    append_meta_and_convert_to_h5ad(extracted_h5_path, dataset_id, args.output)

def parallel_search(args_batches, dataset_id, args_output, max_workers=4):
    with ThreadPoolExecutor(max_workers=min(max_workers, len(args_batches))) as executor:
        futures = {executor.submit(search_and_extract, path, dataset_id, args_output): path for path in args_batches}
        for future in as_completed(futures):
            result = future.result()
            if result:
                return result
    return None

def search_and_extract(tar_path, dataset_id, args_output):
    try:
        with tarfile.open(tar_path, 'r:gz') as tar:
            for member in tar.getmembers():
                if os.path.basename(member.name) == f"{dataset_id}.h5":
                    tar.extract(member, path=args.output)
                    result_path = os.path.join(args_output, member.name)
                    return result_path
    except Exception as e:
        return None
    return None

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

    # Write to .h5ad
    h5ad_path = os.path.join(args_output, f"DISCO-{sample_id}.h5ad")
    adata.write(h5ad_path)
    
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
    parser.add_argument("--dataset", required=True, help="Dataset/sample ID to extract")
    parser.add_argument("--output", required=True, type=str, help="Output directory")
    parser.add_argument('batches', nargs='+', type=Path, help='Paths to batch tar files')
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)