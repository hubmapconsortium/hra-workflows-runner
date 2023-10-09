#!/bin/bash

REFERENCE_DATA_ID=${1:?"A zenodo reference data id must be provided to download!"}
REFERENCE_DATA_DIR=${2:-"$RAW_DATA_DIR/popv/reference-data"}

mkdir -p $REFERENCE_DATA_DIR
zenodo_get $REFERENCE_DATA_ID -o $REFERENCE_DATA_DIR
