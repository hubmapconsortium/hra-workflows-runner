#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

SUMMARIES=$OUTPUT_DIR/sc-transcriptomics-cell-summaries.jsonl.gz
DB=$(realpath $OUTPUT_DIR/dcta-db.duckdb)
SQL=$(realpath $SRC_DIR/create-dcta-db.sql)

rm -f $DB

# Read in raw dataset json files
duckdb $DB -no-stdin -c "CREATE TABLE datasets AS SELECT * FROM read_json('${DATA_REPO_DIR}/*/dataset.json', filename=true, maximum_sample_files=100000)"

# Read in cell summaries
duckdb $DB -no-stdin -c "CREATE TABLE cell_summaries_raw AS SELECT * FROM read_json('${SUMMARIES}', compression='gzip', format='newline_delimited', maximum_object_size=1073741824)"

# Create derived tables
duckdb $DB -no-stdin -init $SQL
