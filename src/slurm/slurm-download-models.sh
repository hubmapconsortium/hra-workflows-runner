#!/bin/bash

#SBATCH -J hra-workflows
#SBATCH -p general
#SBATCH -o slurm-output/download-models/hra-run_%j.txt
#SBATCH -e slurm-output/download-models/hra-run_%j.err
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH -A r00355

module load nodejs
module load python
module load singularity

#Run your program
srun cwl-runner --singularity --no-doc-cache --tmpdir-prefix $TEMP https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/download-models.cwl --outputDirectory $MODELS_DIR
