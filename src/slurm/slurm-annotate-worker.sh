#!/bin/bash

#SBATCH -J hra-workflows-annotate-worker
#SBATCH -A r00355
#SBATCH -p general
#SBATCH -o slurm-output/annotate-worker/hra-run_%j.txt
#SBATCH -e slurm-output/annotate-worker/hra-run_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=4
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mem=32G

module load python
module load singularity

# ---------------------------------------
# Constants
# ---------------------------------------
readonly ALGORITHMS=(azimuth celltypist popv)
readonly CWL_PIPELINE="https://cdn.jsdelivr.net/gh/hubmapconsortium/hra-workflows@main/pipeline.cwl"
readonly DIRS_FILE="$OUTPUT_DIR/annotate-dirs.txt"
readonly MAX_JOBS=4 # Keep in sync with --ntasks
readonly TMP_DIR_PREFIX="$TEMP/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-"

# ---------------------------------------
# Inputs
# ---------------------------------------
declare -i index

readonly FORCE="${FORCE:-"false"}"
readonly RERUN_FAILED="${RERUN_FAILED:-"false"}"
readonly CWL_OPTIONS="${CWL_OPTIONS:-}"
index="${INDEX:-0}"

# ---------------------------------------
# Variables
# ---------------------------------------
declare -a dirs
declare -a locks
declare -a jobs
reschedule="false"

# ----------------------------------------
# Functions
# ----------------------------------------

#######################################
# Cleanup lock files
# Globals:
#   locks
# Arguments:
#   None
#######################################
cleanup() {
  local lock
  for lock in "${locks[@]}"; do
    rm -rf "$lock"
  done
}

#######################################
# Tries to create a lock for a job file
# Globals:
#   locks
# Arguments:
#   Directory of dataset, a path
#   Algorithm, a string
# Returns:
#   0 if the lock was acquired, non-zero otherwise
#######################################
try_lock_job() {
  local -r dir="$1" algorithm="$2"
  local -r file="$dir/job-$algorithm.json.lock"
  local fd

  exec {fd}>"$file"
  if ! flock --nonblock "$fd"; then
    # Make sure to close file descriptor on failure
    exec {fd}>&-
    return 1
  fi

  locks+=("$file")
}

#######################################
# Checks whether a job has run
# Globals:
#   None
# Arguments:
#   Directory of dataset, a path
#   Algorithm, a string
# Outputs:
#   Status of the job if it has run ["success" or "failure"]
# Returns:
#   0 if the job has run, non-zero otherwise
#######################################
check_job_status() {
  local -r dir="$1" algorithm="$2"
  local -r job_dir="$dir/$algorithm"
  local -r report_file="$job_dir/report.json"

  if [[ ! -e "$report_file" ]]; then
    return 1
  fi

  if grep -s success <"$report_file" >/dev/null; then
    echo "success"
  else
    echo "failure"
  fi
}

#######################################
# Checks whether a job should be run based on it's status and global config
# Globals:
#   FORCE
#   RERUN_FAILED
# Arguments:
#   Directory of dataset, a path
#   Algorithm, a string
# Returns:
#   0 if the job should be run, non-zero otherwise
#######################################
check_should_run_job() {
  local status

  if [[ "$FORCE" != "false" ]]; then
    return 0
  fi

  status=$(check_job_status "$1" "$2")
  if [[ "$?" == 0 && ("$status" == "success" || "$RERUN_FAILED" == "false") ]]; then
    return 1
  fi
}

#######################################
# Finds jobs to run
# Globals:
#   ALGORITHMS
#   FORCE
#   MAX_JOBS
#   RERUN_FAILED
#   dirs
#   locks
#   index [Modified]
#   jobs [Modified]
#   reschedule [Modified]
#######################################
find_jobs() {
  local -i count=0 max="$MAX_JOBS"
  local dir algorithm

  for dir in "${dirs[@]}"; do
    for algorithm in "${ALGORITHMS[@]}"; do
      if (("$count" >= "$max")); then
        reschedule="true"
        return 0
      fi

      if ! check_should_run_job "$dir" "$algorithm"; then
        continue
      fi

      if try_lock_job "$dir" "$algorithm"; then
        jobs+=("$dir/job-$algorithm.json")
        count+=1
      fi
    done

    index+=1
  done
}

#######################################
# Run jobs and wait for them to complete
# Globals:
#   CWL_OPTIONS
#   CWL_PIPELINE
#   TMP_DIR_PREFIX
#   jobs
#######################################
run_jobs() {
  local -a pids
  local -i counter=1
  local job dir file

  for job in "${jobs[@]}"; do
    dir=$(dirname "$job")
    file=$(basename "$job")

    pushd "$dir"
    srun cwl-runner --singularity --no-doc-cache --tmpdir-prefix "$TMP_DIR_PREFIX$counter/" "$CWL_OPTIONS" "$CWL_PIPELINE" "$file" &
    pids+=("$!")
    popd

    counter+=1
  done

  wait "${pids[@]}"
}

# ----------------------------------------
# Main logic
# ----------------------------------------

# Print some useful info
printf 'Worker started\n'
printf 'FORCE=%s\n' "$FORCE"
printf 'RERUN_FAILED=%s\n' "$RERUN_FAILED"
printf 'CWL_OPTIONS=%s\n' "$CWL_OPTIONS"
printf 'INDEX=%s\n' "$index"
printf '\n'

# Install cleanup
trap cleanup EXIT

# Read directories
mapfile -t dirs < <(tail -n "+$index" <"$DIRS_FILE")

# Find jobs
find_jobs

# Check if we found any jobs to run
if [[ "${#jobs[@]}" == 0 ]]; then
  printf 'No more jobs. Exiting...\n'
  exit 0
fi

# Run jobs
printf 'Running jobs:\n'
printf '  %s\n' "${jobs[@]}"
run_jobs
printf '\n---------------------\nJobs done!\n---------------------\n'

# Reschedule if there is more work
if [[ "$reschedule" == "true" ]]; then
  # Ensure the latest value of index is exported to the next worker
  readonly sbatch_export="ALL,INDEX=$index"
  sbatch --export "$sbatch_export" "$PROJECT_DIR/src/slurm/slurm-annotate-worker.sh"
fi
