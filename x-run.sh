#!/bin/bash
set -e
source constants.sh

# Link project level sif containers
link_sif .

# Start data download
download_output=$(sbatch $PROJECT_DIR/src/slurm/slurm-run.sh)
download_id=$(get_job_id "$download_output")

if [[ -n "$download_id" ]]; then
  # Start annotations after download
  bash "$PROJECT_DIR/scripts/30x-annotate.sh" --wait-for-job "$download_id" "$@"
fi
