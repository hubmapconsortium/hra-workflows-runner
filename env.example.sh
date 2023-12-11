# Select runner. Available options are:
# "" (local using docker)
# "singularity" (local using singularity)
# "slurm" (slurm workflow manager w/ singularity)
export RUNNER=""

# Must be configured when RUNNER=slurm
# Must be an absolute path
# export PROJECT_DIR="/absolute/path/to/project"

# Highly recommended to be configured when RUNNER=singularity|slurm
# Must be an absolute path
# export SIF_CACHE_DIR="/absolute/path/to/sif/cache"

# Enable if you want to ignore cached files and results
# export FORCE=true

# Enable only specific dataset handlers (by default all handlers are enabled)
# Especially useful when implementing and testing a new handler
# export DATASET_HANDLERS=hubmap,gtex

# Select dataset
export DATASET="sample"
export VERSION="2023-08-15"

# Optionally configure different directory with crosswalking tables
# export CROSSWALKING_TABLES_DIR="./crosswalking/"

# Add handler specific configuration such as tokens, etc.
# See individual handler README.md for available configuration options
export HUBMAP_TOKEN="your token"
export SENNET_TOKEN="your token"
