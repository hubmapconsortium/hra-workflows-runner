import argparse
import json

import os
import sqlite3
import tempfile
import subprocess

import re
import stat
import shutil
import time

from pathlib import Path
import scanpy as sc
import tarfile
import indexed_gzip

import pandas as pd
import math

def main(args):

    # Read metadata
    meta_df = pd.read_csv(args.metadata, sep="\t", dtype=str)
    
    # Process metadata
    meta_dict = {
        row['sample_id']: {k: v for k, v in {**row.to_dict(), 'sample_id': row['sample_id']}.items() 
                          if not (isinstance(v, float) and math.isnan(v))}
        for _, row in meta_df.iterrows()
    }
    # Create the tar map ; Process tar files
    tar_dict = process_tar_files(args.batches)

    # Build sample lists
    tar_list = {k.split('/')[-1].replace('.h5', '') for k in tar_dict}
    meta_list = set(meta_dict.keys())

    # Identify missing tar files from metadata
    missing_tar_from_meta = list(meta_list - tar_list)
    samples_h5ad = [
        sid for sid in missing_tar_from_meta
        if meta_dict.get(sid, {}).get('rds_size') not in [None, '', float('nan')]
    ]
    
    map_dict = build_map_dict(tar_dict, meta_dict, verbose=False)
    samples_list = list(map_dict.keys()) + samples_h5ad

    # print(json.dumps(sorted(samples_list)))

    dataset_id = args.dataset
    sample_id = dataset_id.rstrip('/').split('/')[-1]
    process_sample(sample_id, map_dict, args.batches, args.output)

def process_tar_files(tar_paths):
    tar_dict = {}
    
    for tar_path in tar_paths:                # generic - it is not using for i (1to9) loop
        tar_path = str(tar_path)
        index_path = f"{tar_path}.index.sqlite"
        
        if not Path(tar_path).exists():
            print(f"Skipping missing tar: {tar_path}")
            continue

        # Create unique temporary mount point
        mount_dir = tempfile.mkdtemp(prefix=f"disco_mount_{Path(tar_path).stem}_")

        try:
            # Generate index by mounting the tar file
            result = subprocess.run(
                ['ratarmount', tar_path, mount_dir],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                print(f"Indexing failed for {tar_path}:\n{result.stderr}")
                continue

            # Explicitly unmount
            try:
                if os.name == 'posix':
                    subprocess.run(['umount', '-f', mount_dir], check=True)
                elif os.name == 'nt':
                    subprocess.run(['fusermount', '-u', mount_dir], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Unmount failed: {str(e)}")

        finally:
            safe_cleanup(mount_dir)

        if not Path(index_path).exists():
            print(f"No index created for {tar_path}")
            continue

        # Process the index sqlite files 
        try:
            with sqlite3.connect(index_path) as conn:
                cur = conn.cursor()
                cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = [row[0] for row in cur.fetchall()]
                
                if 'files' not in tables:
                    print(f"'files' table not found in {index_path}")
                    continue

                # Original query and processing logic
                cur.execute("SELECT * FROM files")
                rows = cur.fetchall()
                for row in rows:
                    dir_part = row[0]
                    file_part = row[1]
                    offset = row[2]
                    if dir_part.startswith('/batch_') and file_part.endswith('.h5'):
                        key = f"{dir_part.lstrip('/')}/{file_part}"
                        tar_dict[key] = offset

                # Original per-batch deduplication
                seen = set()
                tar_dict = {k: v for k, v in tar_dict.items() 
                           if (sid := k.split('/')[-1].replace('.h5', '')) not in seen 
                           and not seen.add(sid)}
                
        except sqlite3.Error as e:
            print(f"SQLite error for {tar_path}: {str(e)}")
        finally:
            if Path(index_path).exists():
                Path(index_path).unlink()

    return tar_dict


def build_map_dict(tar_dict, meta_dict, verbose=True):

    map_dict = {}
    tar_sample_map = {}
    for path, offset in tar_dict.items():
        sample_id = path.split('/')[-1].replace('.h5', '')
        tar_sample_map[sample_id] = {
            "path": path,
            "offset": offset
        }

    for sample_id, metadata in meta_dict.items():
        if sample_id in tar_sample_map:
            map_dict[sample_id] = {
                "path": tar_sample_map[sample_id]["path"],
                "offset": tar_sample_map[sample_id]["offset"],
                "metadata": metadata
            }
        else:
            if verbose:
                print(f"Warning: {sample_id} has metadata but no corresponding .h5 file")

    if verbose:
        print(f"Mapped {len(map_dict)} samples")

    return map_dict

def safe_cleanup(path):
    """directory cleanup with permission fixes and retries"""
    for _ in range(3):  # Allow up to 3 retries
        try:
            # Reset permissions recursively
            for root, dirs, files in os.walk(path):
                for name in dirs + files:
                    p = os.path.join(root, name)
                    try:
                        os.chmod(p, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
                    except Exception as e:
                        print(f"Warning: Could not chmod {p}: {str(e)}")
            shutil.rmtree(path)
            return
        except Exception as e:
            print(f"Cleanup failed ({path}), retrying...")
            time.sleep(0.5)
    print(f"Critical: Failed to clean up {path} after 3 attempts")



def extract_h5_from_tar(batch_tar_path, offset, h5_filename, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, h5_filename)
    with indexed_gzip.IndexedGzipFile(str(batch_tar_path)) as f:
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
    tar_path = next(p for p in base_tar_dir if p.name == f"{batch_name}.tar.gz")
    # print("printing tar_path")
    # print(tar_path)
    extracted_h5_path = extract_h5_from_tar(tar_path, offset, h5_filename, output_dir)
    adata = sc.read_10x_h5(extracted_h5_path)

    if metadata:
        for key, value in metadata.items():
            adata.obs[key] = value

    fields = [
        "project_id", "sample_type", "tissue", "anatomical_site", "disease",
        "platform", "age_group", "cell_sorting", "disease_subtype", "treatment",
        "time_point", "subject_id", "age", "gender", "race", "infection", "batch",
        "disease_stage", "genotype", "rna_source", "other_metadata", "disease_grade",
        "cell_number", "median_umi", "rds_md5", "rds_size", "source_cell_line",
        "source_tissue", "source_disease", "source_cell_type", "induced_cell_tissue",
        "collect_time", "sample_id"
    ]

    # Use .iloc[0] to get values for the single sample in adata.obs
    row = adata.obs.iloc[0]

    # Extract only non-empty fields
    filtered_metadata = {
        k: row[k] for k in fields
        if k in row and row[k] not in [None, "", float('nan')]
    }

    # Add cell count and gene count
    filtered_metadata["cell_count"] = len(adata.obs)
    filtered_metadata["gene_count"] = len(adata.var)

    # Print as formatted JSON
    print(json.dumps(filtered_metadata, indent=2))

    h5ad_path = os.path.join(output_dir, f"DISCO-{sample_id}.h5ad")
    adata.write(h5ad_path)


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
