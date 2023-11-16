#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

sif="$SIF_CACHE_DIR/ghcr.io_hubmapconsortium_hra-workflows-runner:main.sif"
work_dir="$PROJECT_DIR"

declare -a bindings
bindings+=("--bind" "$SIF_CACHE_DIR:$SIF_CACHE_DIR")
bindings+=("--bind" "$TEMP:$TEMP")

singularity shell -H "$work_dir" "${bindings[@]}" "$sif"
