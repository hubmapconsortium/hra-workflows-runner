#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

node src/combine-cell-instances.js "$@"
