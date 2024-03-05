#!/bin/bash

# ---------------------------------------
# Constants
# ---------------------------------------
readonly FORCE="${FORCE:-false}"
readonly SKIP_FAILED="${SKIP_FAILED:-false}"

readonly ALGORITHMS=(azimuth celltypist popv)
readonly NUM_ALGORITHMS="${#ALGORITHMS[@]}"

# ---------------------------------------
# Inputs
# ---------------------------------------
readonly CONTROL_FILE="${1:?"Control file not provided"}"
readonly DATASET_DIRS_FILE="${2:?"Dataset directory file not provided"}"
readonly NUM_TASKS="${3:?"Must specify maximum number of runnable jobs to find"}"

# ----------------------------------------
# Functions
# ----------------------------------------

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

# ----------------------------------------
# Main logic
# ----------------------------------------

# Acquire an exclusive lock on the control file
if [[ "$FLOCKER" != "$0" ]]; then
  exec env FLOCKER="$0" flock -e "$CONTROL_FILE" "$0" "$@"
  exit $?
fi

# Read current index
declare -i index
index=$(cat "$CONTROL_FILE")

# Negative index is used to signal that the workers should stop
if ((index < 0)); then
  exit 0
fi

# Read directories
declare -a dirs
declare -i length
mapfile -t dirs <"$DATASET_DIRS_FILE"
length=$((NUM_ALGORITHMS * ${#dirs[@]}))

# Find jobs that have not run or should be rerun
declare -i count=0
declare dir algorithm job_file
while ((index < length && count < NUM_TASKS)); do
  dir="${dirs[$((index / NUM_ALGORITHMS))]}"
  algorithm="${ALGORITHMS[$((index % NUM_ALGORITHMS))]}"

  if should_run "$dir" "$algorithm"; then
    job_file="$dir/job-$algorithm.json"
    if [[ -e "$job_file" ]]; then
      echo "$job_file"
      ((count += 1))
    fi
  fi

  ((index += 1))
done

# Write new index to control file
echo "$index" >"$CONTROL_FILE"
