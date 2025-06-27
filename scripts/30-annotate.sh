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

#######################################
# Checks whether a job has run
# Globals:
#   FORCE
#   SKIP_FAILED
# Arguments:
#   Directory of dataset, a path
#   Algorithm, a string
# Outputs:
#   None
# Returns:
#   0 if the job should run, non-zero otherwise
#######################################
function should_run() {
  local -r report_file="$1/$2/report.json"

  if [[ -e "$report_file" ]]; then
    local -r is_success=$(grep -oe '"status":\s*"success"' "$report_file")
    if [[ -n "$is_success" && "$FORCE" != true ]]; then
      return 1
    elif [[ -z "$is_success" && "$SKIP_FAILED" == true ]]; then
      return 1
    fi
  fi

  return 0
}

# Main logic
if [[ $RUNNER != "slurm" ]]; then
  rm -f jobs.txt jobs2.txt

  for DIR in ${DATASET_DIRS[@]}; do
    for ALGORITHM in azimuth celltypist popv frmatch; do #pan-human-azimuth
      if should_run $DIR $ALGORITHM; then
        if [ -e "${DIR}/job-${ALGORITHM}.json" ]; then
          echo "${SRC_DIR}/run-job.sh ${DIR} ${ALGORITHM}" >> jobs.txt
        fi
      fi
    done
  done

  shuf jobs.txt -o jobs2.txt
  node src/parallel-jobs.js jobs2.txt
  rm -f jobs.txt jobs2.txt
else
  DIRS_FILE="$OUTPUT_DIR/annotate-dirs.txt"
  printf "%s\n" "${DATASET_DIRS[@]}" >$DIRS_FILE

  echo "Use 30x-annotate.sh to run annotations. Exiting..."
  exit $STOP_CODE
fi
