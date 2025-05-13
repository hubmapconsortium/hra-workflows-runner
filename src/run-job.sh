#!/bin/bash
if [ -e .venv ]; then
  source .venv/bin/activate
fi
source constants.sh
shopt -s extglob

DIR=$1
ALGORITHM=$2

CWL_PIPELINE="https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl"
CWL_OPTS=()

if [[ $RUNNER == "slurm" || $RUNNER == "singularity" ]]; then
  CWL_OPTS+=(--singularity)
fi

if [[ -n $TEMP ]]; then
  CWL_OPTS+=(--tmpdir-prefix $TEMP)
fi

cd $DIR
cwl-runner ${CWL_OPTS[@]} ${CWL_PIPELINE} "job-${ALGORITHM}.json" &

wait

exit 0
