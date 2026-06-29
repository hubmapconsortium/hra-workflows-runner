-- CREATE TABLE
--   cell_summaries_raw (
--     "@type" VARCHAR,
--     annotation_method VARCHAR,
--     modality VARCHAR,
--     cell_source VARCHAR,
--     summary STRUCT (cell_id VARCHAR, cell_label VARCHAR, gene_expr JSON, COUNT BIGINT, "@type" VARCHAR, percentage DOUBLE) []
--   );

-- DROP TABLE IF EXISTS cell_summaries;

CREATE TABLE cell_summaries AS
SELECT
  cell_source,
  annotation_method,
  modality,
  summary.cell_id,
  summary.cell_label,
  summary.count,
  summary.percentage,
  json_transform(
    UNNEST(summary.gene_expr ->> '$[*]'),
    '{"gene_id":"VARCHAR","gene_label":"VARCHAR","ensembl_id":"VARCHAR","mean_gene_expr_value":"DOUBLE"}'
  ) AS gene_expr,
  json_transform(
    UNNEST(summary.nsforest_gene_expr ->> '$[*]'),
    '{"gene_id":"VARCHAR","gene_label":"VARCHAR","ensembl_id":"VARCHAR","mean_gene_expr_value":"DOUBLE"}'
  ) AS nsforest_gene_expr
FROM
  (
    SELECT
      cell_source,
      annotation_method,
      modality,
      UNNEST(summary) AS summary
    FROM
      cell_summaries_raw
  );
