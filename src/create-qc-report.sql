CREATE OR REPLACE TABLE qc AS
WITH
    qc_data AS (
        SELECT
            *,
            regexp_extract (filename, '([^/]+)/qc/qc_results/.*', 1) AS folder_name
        FROM
            read_json ('*/qc/qc_results/qc_summary.json', union_by_name = TRUE, filename = TRUE, maximum_sample_files = 100000)
    ),
    datasets AS (
        SELECT
            regexp_extract (filename, '([^/]+)/dataset.json', 1) AS folder_name,
            dataset_id,
            COLUMNS (lambda c: c NOT IN ('dataset_id', 'filename', 'assets', 'rui_location'))
        FROM
            read_json ('*/dataset.json', union_by_name = TRUE, filename = TRUE, maximum_sample_files = 100000)
    )
SELECT
    d.dataset_id,
    q.* EXCLUDE (input_file, filename, thresholds, files),
    d.* EXCLUDE (folder_name, dataset_id)
FROM
    qc_data q
    LEFT JOIN datasets d USING (folder_name);

COPY qc TO 'qc-report.csv.gz';
