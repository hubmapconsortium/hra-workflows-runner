#!/bin/bash
shopt -s extglob
set -e

link_sif() {
  DIR=${1:?"Missing DIR argument to link_sif"}
  if [[ -d $SIF_CACHE_DIR ]]; then
    mkdir -p $DIR

    for SIF in $SIF_CACHE_DIR/*.sif; do
      FILE=$(basename $SIF)
      if [[ ! -e "$DIR/$FILE" ]]; then
        ln -sf $SIF "$DIR/$FILE"
      fi
    done
  fi
}

export -f link_sif
