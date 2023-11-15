#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

if [[ ! -d "$MODELS_DIR" || -z $(ls -A "$MODELS_DIR") ]]; then
  echo "Annotation models are not downloaded. Run '05x-download-models.sh' to download"
  exit 1
fi

# Ensure temp directory exists otherwise mktemp might fail
tmp_dir="${TEMP%'/'}/annotate/"
mkdir -p "$tmp_dir"

start_date=$(date +%Y-%m-%d-%H-%M-%S)
printf 'Annotation run started at %s\n' "$start_date"

control_file=$(mktemp -p "$tmp_dir" "index-$start_date-XXXX.lock")
id="${control_file##*/*-}"
id="${id%.lock}"
printf 'Use "30x-cancel.sh %s" to cancel the run\n' "$id"

echo "0" >"$control_file"
sbatch --array "0-${NUM_WORKERS:-10}" "$SRC_DIR/slurm/slurm-annotate-worker.sh" "$control_file" "$@"
