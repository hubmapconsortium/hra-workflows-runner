#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

MIRROR=https://cdn.humanatlas.io/hra-kg--staging/

for ctann in azimuth celltypist popv frmatch pan-human-azimuth; do
  curl -o crosswalking-tables/${ctann}.csv "${MIRROR}ctann/${ctann}/latest/assets/${ctann}-crosswalk.csv"
done

node src/crosswalks-to-jsonld.js
