#!/bin/bash

#SBATCH -J hra-workflows-annotate-worker
#SBATCH -A r00355
#SBATCH -p general
#SBATCH -o slurm-output/annotate-worker/hra-run_%j.txt
#SBATCH -e slurm-output/annotate-worker/hra-run_%j.err
#SBATCH --mail-type=FAIL,TIME_LIMIT,INVALID_DEPEND,ARRAY_TASKS
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --kill-on-invalid-dep=yes

module load python/3.10.5
module load singularity

# ---------------------------------------
# Constants
# ---------------------------------------
readonly NUM_TASKS="${SLURM_NTASKS:-1}"
readonly HORIZONTAL_LINE='------------------------------------------------'

readonly CWL_PIPELINE_URL="https://cdn.jsdelivr.net/gh/hubmapconsortium/hra-workflows@main/pipeline.cwl"
readonly CWL_TMP_DIR="${TEMP%'/'}/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID"
readonly CWL_OPTIONS="${CWL_OPTIONS:-}"

# ---------------------------------------
# Environment exports
# ---------------------------------------
export APPTAINER_TMPDIR="$TEMP"

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
# Updates the number of active workers
# Globals:
#   SRC_DIR
#   CONTROL_FILE
# Arguments:
#   Number of active workers to add/remove, an integer
#######################################
update_active_workers() {
  bash "$SRC_DIR/slurm/util/update-active-workers.sh" "$CONTROL_FILE" "$1"
}

#######################################
# Finds available jobs
# Globals:
#   SRC_DIR
#   NUM_TASKS
#   CONTROL_FILE
#######################################
find_jobs() {
  bash "$SRC_DIR/slurm/util/find-jobs.sh" "$CONTROL_FILE" "$DATASET_DIRS_FILE" "$NUM_TASKS"
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
  local is_success="" is_ignored_error=""
  
  if [[ -e "$report_file" ]]; then
    is_success="$(grep -oe '"status":\s*"success"' "$report_file")"
    is_ignored_error="$(grep -o -e 'is not supported' -e 'too few unique cells' "$report_file")"
  fi
  if [[ -z "$is_success" && -z "$is_ignored_error" ]]; then
    bash "$SRC_DIR/slurm/util/update-failed-jobs.sh" "$FAILED_DIRS_FILE" "$job_dir"
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
  local tmp_dir

  # Remove the job file path and index
  shift 2

  # Create temp dir
  tmp_dir="$(mktemp -d -p "$CWL_TMP_DIR" "$index-XXXX")"

  # Make a temp dir for apptainer
  export APPTAINER_TMPDIR="${tmp_dir}-apptainer/"
  mkdir -p "${APPTAINER_TMPDIR}"

  # Build cwl-runner command arguments
  args+=(--singularity --no-doc-cache)
  args+=(--tmpdir-prefix "$tmp_dir/") # Note: prefix must end with a '/'
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
declare -i counter=0 worker_count
declare job

# Install traps
trap cleanup EXIT TERM

# Remove worker arguments
shift 3

# Create top temp dir
mkdir -p "$CWL_TMP_DIR"

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
  (run_job "$job" "$counter" "$@") &
  pids+=("$!")
  ((counter += 1))
done

if [[ "${#pids[@]}" == 0 ]]; then
  worker_count="$(update_active_workers -1)"
  if (( "$worker_count" == 0 )); then
    printf 'All workers done. Starting post annotation script\n'
    sbatch "$SRC_DIR/slurm/slurm-post-annotate-run.sh"
  fi

  printf 'No more jobs. Exiting...\n'
  exit 0
fi

# Wait for jobs to finish
wait "${pids[@]}"

# Check if wait terminated due to a signal
if (("$?" > 128)); then
  printf 'Terminated by signal!'
  update_active_workers -1
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
