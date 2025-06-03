import argparse
import os
import sqlite3
import tempfile
import subprocess
import pandas as pd
import json
from pathlib import Path
import re
import stat
import shutil
import time
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

    # Process tar files
    tar_dict = process_tar_files(args.batches)

    # Build sample lists
    tar_list = {k.split('/')[-1].replace('.h5', '') for k in tar_dict}
    meta_list = set(meta_dict.keys())
    
    missing_tar_from_meta = list(meta_list - tar_list)
    samples_h5ad = [
        sid for sid in missing_tar_from_meta
        if meta_dict.get(sid, {}).get('rds_size') not in [None, '', float('nan')]
    ]
    
    map_dict = build_map_dict(tar_dict, meta_dict, verbose=False)
    samples_list = list(map_dict.keys()) + samples_h5ad

    print(json.dumps(sorted(samples_list)))

def process_tar_files(tar_paths):
    tar_dict = {}
    
    for tar_path in tar_paths:
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

        # Process the index sqlite files (original logic)
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


def _get_arg_parser():
    parser = argparse.ArgumentParser(description='Get DISCO dataset listing')
    parser.add_argument('metadata', type=Path, help='Path to metadata TSV file')
    parser.add_argument('batches', nargs='+', type=Path, help='Paths to batch tar files')
    return parser

if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
