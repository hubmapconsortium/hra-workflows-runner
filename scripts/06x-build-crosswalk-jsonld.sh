#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

node src/crosswalks-to-jsonld.js
