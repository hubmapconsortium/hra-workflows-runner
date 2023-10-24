#!/bin/bash

#SBATCH -J hra-workflows
#SBATCH -p general
#SBATCH -o slurm-output/annotate/%A/hra-run_%j_%a.txt
#SBATCH -e slurm-output/annotate/%A/hra-run_%j_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=axbolin@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH -A r00355

module load python
module load singularity

DATASET_DIRS_ARRAY=($DATASET_DIRS)
DIR=${DATASET_DIRS_ARRAY[$SLURM_ARRAY_TASK_ID]}

#Run your program
pushd $DIR
srun cwl-runner --debug --singularity --tmpdir-prefix $TEMP https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl job.json
popd
