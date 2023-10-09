#!/bin/bash

#SBATCH -J hra-workflows
#SBATCH -p general
#SBATCH -o hra-run_%j.txt
#SBATCH -e hra-run_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=axbolin@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH -A r00355

module load python
module load singularity

#Run your program
srun singularity exec -H $PROJECT_DIR --bind $SIF_CACHE_DIR:$SIF_CACHE_DIR $PROJECT_DIR/ghcr.io_hubmapconsortium_hra-workflows-runner\:main.sif ./run.sh
