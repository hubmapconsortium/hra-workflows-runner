#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

for dir in DATA_REPO_DIR/*/; do
  pushd $dir
  cwl-runner --singularity https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl job.yaml
  popd
done
