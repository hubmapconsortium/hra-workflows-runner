#!/usr/bin/env python3

import argparse
import json
import os
import tarfile
import io
import h5py
import anndata
import numpy as np
import sys

ID_COLUMN = os.environ.get("GTEX_COLUMN_ID", "Sample ID")
ORGAN_COLUMN = os.environ.get("GTEX_COLUMN_ORGAN", "Tissue")
SEX_COLUMN = os.environ.get("GTEX_COLUMN_SEX", "Sex")
AGE_COLUMN = os.environ.get("GTEX_COLUMN_AGE", "Age_bin")
DONOR_ID_COLUMN = os.environ.get("GTEX_COLUMN_DONOR", "Participant ID")
TISSUE_SITE_COLUMN = os.environ.get("GTEX_COLUMN_TISSUE_SITE", "Tissue Site Detail")


def find_dataset_location(mapping_file, dataset_id):
    """
    Find which tar file contains the dataset and at what offset.
    
    Args:
        mapping_file (str): Path to the mapping JSON file
        dataset_id (str): ID of the dataset to find
        
    Returns:
        tuple: (tar_file, byte_offset, byte_length) or None if not found
    """
    try:
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            
        # Remove 'disco_' prefix if present in dataset_id
        search_id = dataset_id[6:] if dataset_id.startswith('disco_') else dataset_id
            
        for tar_file, file_data in mapping_data.items():
            for dataset in file_data['datasets']:
                if dataset['name'] == search_id:
                    return (tar_file, dataset['offset'], dataset['length'])
        
        return None
    except Exception as e:
        print(f"Error reading mapping file: {str(e)}", file=sys.stderr)
        sys.exit(1)


def extract_dataset_from_tar(tar_path, byte_offset, byte_length, output_path):
    """
    Extract a specific dataset from a tar.gz file using byte offsets.
    
    Args:
        tar_path (str): Path to the tar.gz file
        byte_offset (int): Starting byte offset in the tar.gz file
        byte_length (int): Length of the dataset in bytes
        output_path (str): Path where to save the extracted dataset
    """
    try:
        with open(tar_path, 'rb') as f:
            # Seek to the byte offset
            f.seek(byte_offset)
            # Read the specified number of bytes
            data = f.read(byte_length)
            
            # Create a file-like object from the bytes
            fileobj = io.BytesIO(data)
            
            # Open as a tar file
            with tarfile.open(fileobj=fileobj, mode='r:gz') as tar:
                # Get the first member (assuming one h5 file per dataset)
                member = tar.next()
                if member is None:
                    raise ValueError("No files found in the tar archive")
                
                # Extract the h5 file to a temporary buffer
                h5_data = tar.extractfile(member).read()
                
                # Write to the output file
                with open(output_path, 'wb') as out_file:
                    out_file.write(h5_data)
                
                # Read the h5 file to get metadata
                with h5py.File(output_path, 'r') as h5f:
                    metadata = {
                        'cell_count': len(h5f['obs_names']) if 'obs_names' in h5f else 0,
                        'gene_count': len(h5f['var_names']) if 'var_names' in h5f else 0
                    }
                    
                    # Try to extract additional metadata if available
                    if 'obs' in h5f:
                        obs = h5f['obs']
                        metadata.update({
                            'tissue': obs['tissue'][()] if 'tissue' in obs else '',
                            'donor_id': obs['donor'][()] if 'donor' in obs else '',
                            'sex': obs['sex'][()] if 'sex' in obs else '',
                            'age': obs['age'][()] if 'age' in obs else ''
                        })
                
                # Print metadata for the caller
                print(f"cell_count: {metadata['cell_count']}")
                print(f"gene_count: {metadata['gene_count']}")
                if 'tissue' in metadata:
                    print(f"tissue: {metadata['tissue']}")
                if 'donor_id' in metadata:
                    print(f"donor_id: {metadata['donor_id']}")
                if 'sex' in metadata:
                    print(f"sex: {metadata['sex']}")
                if 'age' in metadata:
                    print(f"age: {metadata['age']}")
                
    except Exception as e:
        print(f"Error extracting dataset: {str(e)}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Extract a dataset from a tar.gz file")
    parser.add_argument("--dataset", required=True, help="Dataset id")
    parser.add_argument("--mapping-file", required=True, help="Path to the mapping JSON file")
    parser.add_argument("--data-dir", required=True, help="Directory containing the tar.gz files")
    parser.add_argument("--output", required=True, help="Output file path")
    
    args = parser.parse_args()
    
    # Find which tar file contains the dataset
    location = find_dataset_location(args.mapping_file, args.dataset)
    if location is None:
        print(f"Dataset {args.dataset} not found in mapping file", file=sys.stderr)
        sys.exit(1)
        
    tar_file, byte_offset, byte_length = location
    tar_path = os.path.join(args.data_dir, tar_file)
    
    if not os.path.exists(tar_path):
        print(f"Tar file not found: {tar_path}", file=sys.stderr)
        sys.exit(1)
        
    extract_dataset_from_tar(tar_path, byte_offset, byte_length, args.output)


if __name__ == "__main__":
    main()
