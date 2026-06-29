#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

QC_ZIP=$(realpath $OUTPUT_DIR/qc-data.zip)
REPORT=$(realpath $OUTPUT_DIR/qc-report.csv.gz)
REPORT_SQL=$(realpath $SRC_DIR/create-qc-report.sql)

DATASET_DIRS=
for DIR in $(node $SRC_DIR/list-downloaded-dirs.js); do
  DATASET_DIRS+=($DIR)
done

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
    local -r not_supported=$(grep -oE "\"cause\": \"ValueError\('Insufficient cells " "$report_file")
    if [[ ( -n "$is_success" || -n "$not_supported" ) && "$FORCE" != true ]]; then
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
  touch jobs.txt

  for DIR in ${DATASET_DIRS[@]}; do
    for ALGORITHM in qc; do
      if should_run $DIR $ALGORITHM; then
        if [ -e "${DIR}/job-${ALGORITHM}.json" ]; then
          if [ "${MAX_PROCESSES}" == "1" ]; then
            ${SRC_DIR}/run-qc-job.sh ${DIR} ${ALGORITHM}
          else
            echo "${SRC_DIR}/run-qc-job.sh ${DIR} ${ALGORITHM}" >> jobs.txt
          fi
        fi
      fi
    done
  done

  if [ "${MAX_PROCESSES}" != "1" ]; then
    shuf jobs.txt -o jobs2.txt
    node src/parallel-jobs.js jobs2.txt
    rm -f jobs.txt jobs2.txt
  fi
else
  DIRS_FILE="$OUTPUT_DIR/annotate-dirs.txt"
  printf "%s\n" "${DATASET_DIRS[@]}" >$DIRS_FILE

  echo "Use 25-qc-filter.sh to run annotations. Exiting..."
  exit $STOP_CODE
fi


cd $DATA_REPO_DIR

# Create QC report
duckdb -no-stdin -init $REPORT_SQL
mv qc-report.csv.gz $REPORT

# Zip up all qc_results + metadata
find . \( -name 'dataset.json' -o -path '*/qc/qc_results' -o -path '*/qc_results/*' \) -print | zip -@ $QC_ZIP
