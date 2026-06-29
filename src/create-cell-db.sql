-- Create an index for quickly finding unmapped cells
CREATE INDEX unmapped_idx ON cell_instances ((cell_id LIKE 'ASCTB-TEMP:%'));

-- Roll up instances to to cell types
CREATE TABLE cell_types AS
SELECT dataset, organ, tool, cell_id, cell_label, COUNT(*) as cell_count
FROM cell_instances
GROUP BY dataset, organ, tool, cell_id, cell_label
ORDER BY organ, dataset, tool, cell_id;

-- Create an index for quickly finding unmapped cell types
CREATE INDEX unmapped_types_idx ON cell_types ((cell_id LIKE 'ASCTB-TEMP:%'));

-- Create a rolled up unmapped cell types table
CREATE TABLE unmapped_cell_types AS
SELECT *
FROM cell_types
WHERE cell_id LIKE 'ASCTB-TEMP:%';

COPY (SELECT tool, cell_label, cell_id, SUM(cell_count) as count FROM unmapped_cell_types GROUP BY tool, cell_label, cell_id ORDER BY tool, count DESC) TO 'unmapped_cell_types.csv';
