import argparse
import os
import json
import re
from pathlib import Path
import scanpy as sc
import tarfile
import indexed_gzip


def extract_h5_from_tar(batch_tar_path, offset, h5_filename, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, h5_filename)
    with indexed_gzip.IndexedGzipFile(batch_tar_path) as f:
        f.seek(offset)
        with tarfile.open(fileobj=f, mode='r|') as tar:
            for member in tar:
                if member.name.endswith(h5_filename):
                    tar.extract(member, path=output_dir)
                    extracted_path = os.path.join(output_dir, member.name)
                    if extracted_path != output_path:
                        os.rename(extracted_path, output_path)
                    return output_path
    raise FileNotFoundError(f"{h5_filename} not found in {batch_tar_path} at offset {offset}")


def process_sample(sample_id, disco_map, base_tar_dir, output_dir):
    entry = disco_map.get(sample_id)
    if not entry:
        raise ValueError(f"{sample_id} not found in disco_map")

    path = entry['path']
    offset = entry['offset']
    metadata = entry['metadata']

    batch_name = path.split('/')[0]
    h5_filename = path.split('/')[-1]
    tar_path = os.path.join(base_tar_dir, f"{batch_name}.tar.gz")

    extracted_h5_path = extract_h5_from_tar(tar_path, offset, h5_filename, output_dir)
    adata = sc.read_10x_h5(extracted_h5_path)

    if metadata:
        for key, value in metadata.items():
            adata.obs[key] = value

    print(json.dumps({
        "organ": metadata.get("tissue", ""),
        "sex": metadata.get("gender", ""),
        "age": metadata.get("age", ""),
        "donor_id": metadata.get("subject_id", ""),
        "cell_count": len(adata.obs),
        "gene_count": len(adata.var),
        "tissue_site": metadata.get("anatomical_site", "")
    }, indent=2))

    h5ad_path = os.path.join(output_dir, f"{sample_id}.h5ad")
    adata.write(h5ad_path)


def main(args):
    with open(args.map_file, 'r') as f:
        disco_map = json.load(f)
    dataset_id = args.dataset
    sample_id = dataset_id.rstrip('/').split('/')[-1]
    process_sample(sample_id, disco_map, args.tar_dir, args.output)


def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Extract DISCO dataset from .tar.gz batches")
    parser.add_argument("--dataset", required=True, help="Dataset/sample ID to extract")
    parser.add_argument("--output", required=True, type=str, help="Output directory")
    parser.add_argument("--map_file", required=True, type=Path, help="Path to disco_map.json")
    parser.add_argument("--tar_dir", required=True, type=Path, help="Directory containing batch .tar.gz files")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
