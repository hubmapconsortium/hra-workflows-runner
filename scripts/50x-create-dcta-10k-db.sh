#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

DB=$(realpath $OUTPUT_DIR/dcta-db-10k.duckdb)
SQL=$(realpath $SRC_DIR/create-dcta-db.sql)

rm -f $DB

# Read in cell summaries
duckdb $DB -no-stdin -c "CREATE TABLE cell_summaries_raw AS SELECT unnest(\"@graph\"[1]) FROM read_json('${DATA_REPO_DIR}/*/*/summary.jsonld', ignore_errors=true, auto_detect=true, union_by_name=true, maximum_object_size=1073741824)"

# Read in raw dataset json files
duckdb $DB -no-stdin -c "CREATE TABLE datasets AS SELECT * FROM read_json('${DATA_REPO_DIR}/*/dataset.json', filename=true, maximum_sample_files=100000)"

# Create derived tables
duckdb $DB -no-stdin -init $SQL
