#!/bin/bash

#SBATCH -J hra-workflows
#SBATCH -p general
#SBATCH -o slurm-output/run/hra-run_%j.txt
#SBATCH -e slurm-output/run/hra-run_%j.err
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH -A r00355

module load python
module load singularity

#Run your program
srun singularity exec -H $PROJECT_DIR --bind $SIF_CACHE_DIR:$SIF_CACHE_DIR $PROJECT_DIR/ghcr.io_hubmapconsortium_hra-workflows-runner\:main.sif ./run.sh
