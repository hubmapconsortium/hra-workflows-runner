shopt -s expand_aliases

# Set global configuration
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export PYTHONPATH=".:./src"
export GPG_TTY=$(tty)

# Load environment
source env.sh

# Shorthands and configuration options
export DATASETS_DIR=${DATASETS_DIR:="./datasets"}
export RAW_DATA_DIR=${RAW_DATA_DIR:="./raw-data"}
