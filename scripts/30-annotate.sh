#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CWL_PIPELINE="https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl"
CWL_OPTS=()
readarray -t DATASET_DIRS < <(node $SRC_DIR/list-downloaded-dirs.js)

if [[ $RUNNER == "slurm" || $RUNNER == "singularity" ]]; then
  CWL_OPTS+=(--singularity)
fi

if [[ -n $TEMP ]]; then
  CWL_OPTS+=(--tmpdir-prefix $TEMP)
fi

if [[ ${CWL_OPTS[@]} =~ "singularity" ]]; then
  for DIR in ${DATASET_DIRS[@]}; do
    link_sif $DIR
  done
fi

# Main logic
if [[ $RUNNER != "slurm" ]]; then
  for DIR in ${DATASET_DIRS[@]}; do
    pushd $DIR
    cwl-runner ${CWL_OPTS[@]} $CWL_PIPELINE job.json
    popd
  done
else
  SBATCH_FILE="$OUTPUT_DIR/sbatch.sh"
  EXPORTS="ALL,DATASET_DIRS=\"${DATASET_DIRS[@]}\""
  RANGE="0-$((${#DATASET_DIRS[@]} - 1))"

  echo "sbatch --array $RANGE --export $EXPORTS $SRC_DIR/slurm/slurm-annotate.sh" >$SBATCH_FILE
  chmod +x $SBATCH_FILE

  echo "Generated sbatch commands. Use 30x-annotate.sh to run annotations. Exiting..."
  exit 1
fi
