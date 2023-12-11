#!/bin/bash
set -e
source constants.sh

link_sif .
sbatch $PROJECT_DIR/src/slurm/slurm-run.sh
