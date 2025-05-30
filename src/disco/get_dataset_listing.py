#!/usr/bin/env python3

import json
import sys
import os

def get_dataset_listing(mapping_file):
    """
    Read the mapping file that contains information about which datasets are in which tar.gz files
    and at what byte offsets.
    
    Returns a list of dataset information in JSON format.
    """
    try:
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            
        datasets = []
        for tar_file, file_data in mapping_data.items():
            for dataset in file_data['datasets']:
                dataset_info = {
                    'id': f"disco_{dataset['name']}",
                    'tar_file': tar_file,
                    'byte_offset': dataset['offset'],
                    'byte_length': dataset['length'],
                    'tissue': dataset.get('tissue', ''),
                    'donor_id': dataset.get('donor_id', ''),
                    'cell_count': dataset.get('cell_count', 0),
                }
                datasets.append(dataset_info)
                
        return datasets
    except Exception as e:
        print(f"Error reading mapping file: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 get_dataset_listing.py <mapping_file>", file=sys.stderr)
        sys.exit(1)
        
    mapping_file = sys.argv[1]
    if not os.path.exists(mapping_file):
        print(f"Mapping file not found: {mapping_file}", file=sys.stderr)
        sys.exit(1)
        
    datasets = get_dataset_listing(mapping_file)
    print(json.dumps(datasets))

if __name__ == "__main__":
    main()
