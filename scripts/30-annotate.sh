#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

readarray -t DATASET_DIRS < <(node $SRC_DIR/list-downloaded-dirs.js)

if [[ $RUNNER == "slurm" || $RUNNER == "singularity" ]]; then
  for DIR in ${DATASET_DIRS[@]}; do
    link_sif $DIR
  done
fi

# Main logic
if [[ $RUNNER != "slurm" ]]; then
  rm -f progress
  for DIR in ${DATASET_DIRS[@]}; do
    echo $(basename $DIR) >> progress
    for ALGORITHM in azimuth celltypist popv; do
      queue_job "${SRC_DIR}/run-job.sh ${DIR} ${ALGORITHM}"
    done
  done

  wait_for_empty_queue
else
  DIRS_FILE="$OUTPUT_DIR/annotate-dirs.txt"
  printf "%s\n" "${DATASET_DIRS[@]}" >$DIRS_FILE

  echo "Use 30x-annotate.sh to run annotations. Exiting..."
  exit $STOP_CODE
fi
