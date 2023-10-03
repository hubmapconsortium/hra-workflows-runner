#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CWL_PIPELINE="https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl"
CWL_OPTS=()

if [[ $RUNNER == "slurm" || $RUNNER == "singularity" ]]; then
  CWL_OPTS+=(--singularity)

  if [[ -d $SIF_CACHE_DIR ]]; then
    LINK_SIF_CACHE=true
  fi
fi

if [[ -n $TEMP ]]; then
  CWL_OPTS+=(--tmpdir-prefix $TEMP)
fi

link_sif_cache() {
  for SIF in $SIF_CACHE_DIR/*.sif; do
    FILE=$(basename $SIF)
    ln -sf $SIF "$1/$FILE"
  done
}

for DIR in $(node $SRC_DIR/list-downloaded-dirs.js); do
  pushd $DIR

  if [[ -n $LINK_SIF_CACHE ]]; then
    link_sif_cache $DIR
  fi

  if [[ $RUNNER == "slurm" ]]; then
    sbatch $PROJECT_DIR/src/slurm/slurm-annotate.sh
  else
    cwl-runner ${CWL_OPTS[@]} $CWL_PIPELINE job.json
  fi

  popd
done
