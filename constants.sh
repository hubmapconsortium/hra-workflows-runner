shopt -s expand_aliases

# Set global configuration
export LANG=C.UTF-8
export PYTHONPATH=".:./src"
export GPG_TTY=$(tty)

# Load environment
if [ -e env.sh ]; then
  source env.sh
fi

# Check required configuration
_=${DATASET:?"No dataset selected!"}
_=${VERSION:?"No version selected!"}

# Shorthands and configuration options
export RUNNER=${RUNNER:-"docker"}
export DATASETS_DIR=${DATASETS_DIR:-"./datasets"}
export CROSSWALKING_TABLES_DIR=${CROSSWALKING_TABLES_DIR:-"./crosswalking-tables"}
export RAW_DATA_DIR=${RAW_DATA_DIR:-"./raw-data"}
export OUTPUT_DIR="$RAW_DATA_DIR/$DATASET/$VERSION"
export DATA_REPO_DIR="$RAW_DATA_DIR/data-repo"
export CACHE_DIR="$RAW_DATA_DIR/tmp"
export MODELS_DIR="$RAW_DATA_DIR/models"
export SRC_DIR="./src"

# Absolute path to project. Used when starting slurm jobs
export PROJECT_DIR=${PROJECT_DIR:-"/N/project/hra/hra-workflows-runner"}
export SIF_CACHE_DIR=${SIF_CACHE_DIR:-"/N/project/hra/sif-cache"}
export TEMP=${TEMP:-"$PROJECT_DIR/tmp/"}

# Configure slurm
export SLURM_UMASK="ug+rw"

# Load dataset configuration
source $DATASETS_DIR/$DATASET/config.sh

# Check required dataset configuration
_=${DATASET_COLUMN_ID:?"No dataset list identifier column set!"}

# Load shared functions
for SCRIPT in $SRC_DIR/shell_util/*.sh; do
  source $SCRIPT
done
