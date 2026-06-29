#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

INSTANCES=$OUTPUT_DIR/sc-transcriptomics-cell-instances.csv.gz
DB=$(realpath $OUTPUT_DIR/cell-db.duckdb)
SQL=$(realpath $SRC_DIR/create-cell-db.sql)

rm -f $DB

# Read in cell instances
zcat $INSTANCES | csvformat -T | duckdb $DB -c "CREATE TABLE cell_instances AS SELECT * FROM read_csv('/dev/stdin')"

# Create derived tables
cd $OUTPUT_DIR
duckdb $DB -no-stdin -init $SQL
