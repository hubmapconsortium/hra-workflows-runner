#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

node src/download.js
find "$DATA_REPO_DIR" -empty -type d -delete
