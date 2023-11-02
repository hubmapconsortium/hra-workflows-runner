#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

if [[ ! -d "$MODELS_DIR" || -z $(ls -A "$MODELS_DIR") ]]; then
  echo "Annotation models are not downloaded. Run '05x-download-models.sh' to download"
  exit 1
fi

sbatch --array "0-${NUM_WORKERS:-10}" "$SRC_DIR/slurm/slurm-annotate-worker.sh"
