#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

for dir in DATA_REPO_DIR/*/; do
  pushd $dir
  sbatch $SRC_DIR/src/slurm/slurm-annotate.sh
  popd
done
