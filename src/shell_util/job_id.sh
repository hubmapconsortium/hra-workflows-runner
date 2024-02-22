#!/bin/bash
shopt -s extglob
set -e

get_job_id() {
  declare -r sbatch_output="${1:?Missing slurm output}"
  declare -r id_regex="job[[:space:]]+([0-9]+)"

  if [[ "$sbatch_output" =~ $id_regex ]]; then
    echo "${BASH_REMATCH[1]}"
  fi
}

export -f get_job_id
