#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CONTAINERS=(azimuth celltypist popv pan-human-azimuth frmatch extract-summary gene-expression crosswalking nsforest) 

mkdir -p "$SIF_CACHE_DIR"
cd "$SIF_CACHE_DIR"

# Main hra-workflows-runner container
docker pull ghcr.io/hubmapconsortium/hra-workflows-runner:main

# hra-workflows containers
for name in "${CONTAINERS[@]}"; do
  url="ghcr.io/hubmapconsortium/hra-workflows/$name:main"
  docker pull $url
done

# nodejs
docker pull node:20-alpine
