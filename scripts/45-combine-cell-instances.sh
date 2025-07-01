#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

MAX_PROCESSES=`nproc` node src/combine-cell-instances.js "$@"
