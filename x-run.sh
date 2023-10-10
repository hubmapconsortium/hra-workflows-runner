#!/bin/bash
set -e
source constants.sh

sbatch $PROJECT_DIR/src/slurm/slurm-run.sh
