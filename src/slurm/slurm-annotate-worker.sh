#!/bin/bash

#SBATCH -J hra-workflows-annotate-worker
#SBATCH -A r00355
#SBATCH -p general
#SBATCH -o slurm-output/annotate-worker/hra-run_%j.txt
#SBATCH -e slurm-output/annotate-worker/hra-run_%j.err
#SBATCH --mail-type=FAIL,TIME_LIMIT,INVALID_DEPEND,ARRAY_TASKS
#SBATCH --ntasks=4
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mem=128G
#SBATCH --kill-on-invalid-dep=yes

module load python
module load singularity

# ---------------------------------------
# Constants
# ---------------------------------------
readonly NUM_TASKS="${SLURM_NTASKS:-4}"
readonly HORIZONTAL_LINE='------------------------------------------------'

readonly CWL_PIPELINE_URL="https://cdn.jsdelivr.net/gh/hubmapconsortium/hra-workflows@main/pipeline.cwl"
readonly CWL_TMP_DIR="${TEMP%'/'}/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID"
readonly CWL_OPTIONS="${CWL_OPTIONS:-}"

# ---------------------------------------
# Inputs
# ---------------------------------------
readonly CONTROL_FILE="${1:?"Missing control file"}"
readonly DATASET_DIRS_FILE="${2:?"Missing dataset directory file"}"
readonly FAILED_DIRS_FILE="${3:?"Missing failed dataset output file"}"

# ----------------------------------------
# Functions
# ----------------------------------------

#######################################
# Perform cleanup
# Globals:
#   CWL_TMP_DIR
#   pids
#######################################
cleanup() {
  # Stop any annotation tools still running
  if [[ -n "$pids" ]]; then
    kill -s TERM --timeout 10000 KILL "${pids[@]}" &>/dev/null
  fi

  # Clean up after cwl in case it failed to do so properly itself
  rm -rf "$CWL_TMP_DIR"
}

#######################################
# Finds available jobs
# Globals:
#   SRC_DIR
#   NUM_TASKS
#   CONTROL_FILE
#######################################
find_jobs() {
  exec "$SRC_DIR/slurm/util/find-jobs.sh" "$CONTROL_FILE" "$DATASET_DIRS_FILE" "$NUM_TASKS"
}

#######################################
# Checks the result of a job run
# Globals:
#   FAILED_DIRS_FILE
# Arguments:
#   Job file path, a string
#######################################
check_job_result() {
  local -r job="$1"
  local -r job_dir="$(dirname "$job")"

  [[ "$job" =~ job-(.+)\.json ]]
  local -r algorithm="${BASH_REMATCH[1]}"
  local -r report_file="$job_dir/$algorithm/report.json"
  local is_success="" is_unsupported_organ=""
  
  if [[ -e "$report_file" ]]; then
    is_success="$(grep -oe '"status":\s*"success"' "$report_file")"
    is_unsupported_organ="$(grep -oe 'is not supported' "$report_file")"
  fi
  if [[ -z "$is_success" && -z "$is_unsupported_organ" ]]; then
    exec "$SRC_DIR/slurm/util/update-failed-jobs.sh" "$FAILED_DIRS_FILE" "$job_dir"
  fi
}

#######################################
# Checks whether a job has run
# Globals:
#   CWL_PIPELINE_URL
#   CWL_TMP_DIR
#   CWL_OPTIONS
# Arguments:
#   Job file path, a string
#   Job index, an integer
#######################################
run_job() {
  local -r job="$1" index="$2"
  local -a args

  # Remove the job file path and index
  shift 2

  # Build cwl-runner command arguments
  args+=(--singularity --no-doc-cache)
  args+=(--tmpdir-prefix "$CWL_TMP_DIR/$index/")
  args+=("$@")
  args+=("$CWL_PIPELINE_URL")
  args+=("$(basename "$job")")

  # Run command
  pushd "$(dirname "$job")" >/dev/null || exit
  srun time cwl-runner "${args[@]}"
  popd >/dev/null || exit
  check_job_result "$job"
}

#######################################
# Reschedule another worker
# Globals:
#   SRC_DIR
# Arguments:
#   The same arguments that was passed to this script
#######################################
reschedule_worker() {
  printf 'Rescheduling more work\n'
  sbatch "$SRC_DIR/slurm/slurm-annotate-worker.sh" "$@"
}

# ----------------------------------------
# Main logic
# ----------------------------------------

declare -a pids
declare -i counter=0
declare job

# Install traps
trap cleanup EXIT TERM

# Remove worker arguments
shift 3

cat <<EOF
$HORIZONTAL_LINE
  Worker started!
  $(date)
  DATASET: $DATASET - $VERSION
  CONTROL_FILE=$CONTROL_FILE
  ARGS=$@
  FORCE=${FORCE:-false}
  SKIP_FAILED=${SKIP_FAILED:-false}
$HORIZONTAL_LINE
EOF

# Find and start jobs
for job in $(find_jobs); do
  printf 'Starting job: %s\n' "$job"
  run_job "$job" "$counter" "$@" &
  pids+=("$!")
  ((counter += 1))
done

if [[ "${#pids[@]}" == 0 ]]; then
  printf 'No more jobs. Exiting...\n'
  exit 0
fi

# Wait for jobs to finish
wait "${pids[@]}"

# Check if wait terminated due to a signal
if (("$?" > 128)); then
  printf 'Terminated by signal!'
  exit 0
fi

cat <<EOF
$HORIZONTAL_LINE
  Jobs done!
  $(date)
$HORIZONTAL_LINE
EOF

# Reschedule another worker
reschedule_worker "$CONTROL_FILE" "$DATASET_DIRS_FILE" "$FAILED_DIRS_FILE" "$@"
