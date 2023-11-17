#!/bin/bash
source constants.sh
shopt -s extglob
set -e

# ---------------------------------------
# Argument parsing
# ---------------------------------------
declare -a args worker_args
declare -i num_workers="${NUM_WORKERS:-10}"

while (("$#" > 0)); do
  case $1 in
  -n | --num-workers)
    num_workers="$2"
    shift 2
    ;;
  --)
    shift              # Remove --
    worker_args=("$@") # All arguments after -- are passed to the worker script
    shift "$#"         # Consume all arguments
    ;;

  -* | --*)
    echo "Unknown option $1"
    exit 1
    ;;

  *)
    args+=("$1")
    shift
    ;;
  esac
done

set -- "${args[@]}" # Restore positional arguments

if (("$#" > 1)); then
  echo "Unknown extra positional arguments ${@:1}"
  echo "Use '-- arg1 arg2 .. argN' to pass arguments to the worker script"
  exit 1
fi

# ---------------------------------------
# Check preconditions
# ---------------------------------------
if [[ ! -d "$MODELS_DIR" || -z $(ls -A "$MODELS_DIR") ]]; then
  echo "Annotation models are not downloaded. Run '05x-download-models.sh' to download"
  exit 1
fi

# ---------------------------------------
# Setup control file
# ---------------------------------------
declare -a matches
declare id start_date control_file tmp_dir

# Ensure temporary directory exists
tmp_dir="${TEMP%'/'}/annotate/"
mkdir -p "$tmp_dir"

if [[ -z $1 ]]; then
  # Create a new control file
  start_date=$(date +%Y-%m-%d-%H-%M-%S)
  control_file=$(mktemp -p "$tmp_dir" "index-$start_date-XXXX.lock")

  # Initialize control file content
  echo "0" >"$control_file"
else
  # Try to find an existing control file
  mapfile -t matches < <(find "$tmp_dir" -name "index-*$1*.lock")
  control_file="${matches[0]}"

  if (("${#matches[@]}" == 0)); then
    echo "No matching control file for id: '$1'"
    exit 1
  elif (("${#matches[@]}" > 1)); then
    echo "Multiple matching control files: ${matches[@]}"
    exit 1
  fi
fi

id="${control_file##*/*-}"
id="${id%.lock}"

# ---------------------------------------
# Start worker
# ---------------------------------------
echo "Starting workers"
echo "Control file: $control_file"
echo "Number of new workers: $num_workers"
echo "Use '30x-annotate $id' to start additional workers"
echo "Use '30x-cancel.sh $id' to cancel the run"
sbatch --array "1-$num_workers" "$SRC_DIR/slurm/slurm-annotate-worker.sh" "$control_file" "${worker_args[@]}"
