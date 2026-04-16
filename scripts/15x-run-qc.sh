#!/bin/bash
source constants.sh

QC_ZIP=$(realpath $OUTPUT_DIR/qc-data.zip)
REPORT=$(realpath $OUTPUT_DIR/qc-report.csv.gz)
REPORT_SQL=$(realpath $SRC_DIR/create-qc-report.sql)

# Run QC on all datasets
MAX_QC_PASSES=3
for pass in $(seq 1 $MAX_QC_PASSES); do
  echo "QC pass ${pass}..."
  rm -f qc-jobs.txt
  while IFS= read -r data_file; do
    f=${data_file%/data.h5ad}
    if [ ! -e "$f/qc_results/qc_summary.json" ]; then
      echo "./src/run-qc-job.sh $f" >> qc-jobs.txts
    fi
  done < <(find raw-data-v1.1/data-repo -name 'data.h5ad')

  if [ ! -s qc-jobs.txt ]; then
    echo "No QC jobs remaining (pass $pass/$MAX_QC_PASSES)."
    break
  fi

  # Shuffle qc-jobs.txt in-place
  mapfile -t lines < qc-jobs.txt; printf '%s\n' "${lines[@]}" | shuf > qc-jobs.txt

  node src/parallel-jobs.js qc-jobs.txt
done

# Create QC report
duckdb -no-stdin -init $REPORT_SQL
mv qc-report.csv.gz $REPORT

# Zip up all qc_results + metadata
cd raw-data-v1.1
find data-repo \( -name 'dataset.json' -o -path '*/qc_results' -o -path '*/qc_results/*' \) -print | zip -@ $QC_ZIP
