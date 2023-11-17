#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

node src/combine-summary-jsonld.js "$@"
