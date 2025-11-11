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
  # If singularity doesn't exist, do not add the singularity option
  if [[ `which singularity` != "" ]]; then
    CWL_OPTS+=(--singularity)
  fi
fi

if [[ -n $TEMP ]]; then
  if [[ `which singularity` != "" ]]; then
    CWL_OPTS+=(--tmpdir-prefix $TEMP)
  else
    # Docker needs to have temp data in /tmp
    CWL_OPTS+=(--tmpdir-prefix /tmp/dcta/tmp)
  fi
fi

cd $DIR
cwl-runner ${CWL_OPTS[@]} ${CWL_PIPELINE} "job-${ALGORITHM}.json" &

wait

exit 0
