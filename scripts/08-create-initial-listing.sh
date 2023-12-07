#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

if [[ ! -e "$DATASETS_DIR/$DATASET/listing.csv" || "$FORCE" == true ]]; then
  node src/create-initial-listing.js
fi
