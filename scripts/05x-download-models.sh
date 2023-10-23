#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CWL_PIPELINE="https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/download-models.cwl"
CWL_OPTS=()

if [[ $RUNNER == "slurm" || $RUNNER == "singularity" ]]; then
  CWL_OPTS+=(--singularity)
fi

if [[ ${CWL_OPTS[@]} =~ "singularity" ]]; then
  link_sif .
fi

if [[ $RUNNER != "slurm" ]]; then
  cwl-runner ${CWL_OPTS[@]} $CWL_PIPELINE --outputDirectory $MODELS_DIR
else
  sbatch $SRC_DIR/slurm/slurm-download-models.sh
fi
