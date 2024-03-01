#!/bin/bash

# ---------------------------------------
# Inputs
# ---------------------------------------
readonly FAILED_DIRS_FILE="${1:?"Missing failed dataset output file"}"
readonly FAILED_DIR="${2:?"Missing failed dataset directory"}"

# ----------------------------------------
# Main logic
# ----------------------------------------

# Acquire an exclusive lock on the directory file
if [[ "$FLOCKER" != "$0" ]]; then
  exec env FLOCKER="$0" flock -e "$FAILED_DIRS_FILE" "$0" "$@"
  exit $?
fi

# Create if it doesn't exist
if [[ ! -e "$FAILED_DIRS_FILE" ]]; then
  touch "$FAILED_DIRS_FILE"
fi

# Append directory only if not already in the file
exists_in_file="$(grep -F "$FAILED_DIR" "$FAILED_DIRS_FILE")"
if [[ -z "$exists_in_file" ]]; then
  echo "$FAILED_DIR" >>"$FAILED_DIRS_FILE"
fi
