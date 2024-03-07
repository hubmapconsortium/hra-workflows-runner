#!/bin/bash

# ---------------------------------------
# Inputs
# ---------------------------------------
readonly CONTROL_FILE="${1:?"Missing active worker control file"}"
readonly WORKER_CONTROL_FILE="${CONTROL_FILE/"index"/"active-workers"}"
readonly CHANGE_AMOUNT="${2:?"Missing change amount"}"

# ----------------------------------------
# Main logic
# ----------------------------------------

# Acquire an exclusive lock on the control file
if [[ "$FLOCKER" != "$0" ]]; then
  exec env FLOCKER="$0" flock -e "$WORKER_CONTROL_FILE" "$0" "$@"
  exit $?
fi

# Create if it doesn't exist
if [[ ! -e "$WORKER_CONTROL_FILE" ]]; then
  touch "$WORKER_CONTROL_FILE"
fi

# Read current count
declare -i count
count="$(cat "$WORKER_CONTROL_FILE")"
count="${count:-0}"

# Update and output new count both to file and stdout
echo "$((count + CHANGE_AMOUNT))" | tee "$WORKER_CONTROL_FILE"
