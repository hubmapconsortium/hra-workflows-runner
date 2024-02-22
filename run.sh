#!/bin/bash
set -e
source constants.sh

# Code to signal an early but normal exit
export STOP_CODE=99

echo "Run started on $(date)..."
echo
for f in scripts/??-*.sh; do
  echo "Running $f..."
  if time bash $f; then
    echo
  elif [[ "$?" == "$STOP_CODE" ]]; then
    echo
    echo "Early exit signaled from $f"
    exit
  else
    exit "$?"
  fi
done
