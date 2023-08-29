#!/bin/bash

#SBATCH -J hra-workflows
#SBATCH -p general
#SBATCH -o hra-run_%j.txt
#SBATCH -e hra-run_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=axbolin@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH -A general

module load python
module load singularity

#Run your program
srun cwl-runner --singularity https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl job.yaml