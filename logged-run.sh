#!/bin/bash
source constants.sh

LOG=log.txt

time bash -c "time ./run.sh --from 00 --to 99 2>&1" | tee $LOG
