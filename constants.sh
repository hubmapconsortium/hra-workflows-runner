shopt -s expand_aliases

# Set global configuration
export LC_ALL=C.UTF-8
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
export DATASETS_DIR=${DATASETS_DIR:="./datasets"}
export RAW_DATA_DIR=${RAW_DATA_DIR:="./raw-data"}
export OUTPUT_DIR="$RAW_DATA_DIR/$DATASET/$VERSION"
export DATA_REPO_DIR="$RAW_DATA_DIR/data-repo"
export CACHE_DIR="$RAW_DATA_DIR/tmp"
export SRC_DIR="./src"

# Absolute path to project. Used when starting slurm jobs
export PROJECT_DIR="/N/project/hra/hra-workflows-runner"

# Load dataset configuration
source $DATASETS_DIR/$DATASET/config.sh

# Check required dataset configuration
_=${DATASET_LIST_URL:?"No dataset list url set!"}
_=${DATASET_COLUMN_ID:?"No dataset list identifier column set!"}
