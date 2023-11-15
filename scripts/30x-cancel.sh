#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

id="${1:?"Must provide the run id to cancel a run"}"

tmp_dir="${TEMP%'/'}/annotate"
control_file=$(find "$tmp_dir" -name "*$id.lock")

if [[ ! -e "$control_file" ]]; then
  echo "Could not found control file for run $id"
  exit 1
fi

(
  flock 9 || exit 1
  printf 'Approximate number of executed jobs: %s\n' "$(<$control_file)"
  echo "-1" >"$control_file"
) 9<"$control_file"
