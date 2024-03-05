#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CONTAINERS=(azimuth celltypist popv extract-summary gene-expression crosswalking)

mkdir -p "$SIF_CACHE_DIR"
cd "$SIF_CACHE_DIR"

# Main hra-workflows-runner container
sif="ghcr.io_hubmapconsortium_hra-workflows-runner:main.sif"
url="docker://ghcr.io/hubmapconsortium/hra-workflows-runner:main"
apptainer pull --force --name "$sif" "$url"

# hra-workflows containers
for name in "${CONTAINERS[@]}"; do
  sif="ghcr.io_hubmapconsortium_hra-workflows_$name:main.sif"
  url="docker://ghcr.io/hubmapconsortium/hra-workflows/$name:main"
  apptainer pull --force --name "$sif" "$url"
done

# nodejs
sif="node_alpine.sif"
url="docker://node:20-alpine"
apptainer pull --force --name "$sif" "$url"
