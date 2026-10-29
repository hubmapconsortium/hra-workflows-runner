/**
 * Worker thread that reads a cell annotation CSV file for a specific dataset/tool,
 * computes a confidence score for each row, formats the results as CSV lines, and
 * streams each line back to the main thread using `parentPort.postMessage()`.
 *
 * This script is invoked by the main thread with the following `workerData`:
 * {
 *   tool: string,
 *   dataset: object,
 *   popvMethodCount: number,
 *   path: string,
 *   cellIdColumn?: string,
 *   cellLabelColumn?: string
 * }
 */

import { parentPort, workerData } from 'node:worker_threads';
import Papa from 'papaparse';
import { readCsv } from './util/csv.js';

const { tool, dataset, popvMethodCount, path, cellIdColumn, cellLabelColumn } = workerData;

try {
  const csvStream = readCsv(path);

  for await (const cell of csvStream) {
    /**
     * Compute confidence score using available fields, with fallbacks.
     * Divides popv_prediction_score by the method count if necessary.
     */
    const confidence_score =
      parseFloat(cell['mapping.score']) || // azimuth
      parseFloat(cell.conf_score) || // celltypist
      parseFloat(cell.final_level_confidence) || // pan-human-azimuth
      parseFloat(cell.frmatch_confidence) || // frmatch
      parseFloat(cell.popv_prediction_score) / popvMethodCount || // popv
      0;

    /**
     * Format one entry for this cell. Order matches the output header.
     */
    const entry = [
      dataset.dataset_id,
      dataset.organ ?? '',
      tool,
      cell[''] || cell['index'],
      cellIdColumn ? cell[cellIdColumn] : cell.clid,
      cellLabelColumn ? cell[cellLabelColumn] : cell.hra_prediction,
      // Author provided ids are already CL ids
      cellIdColumn ? 'skos:exactMatch' : cell.match_type,
      confidence_score,
    ];

    // Emit one CSV line back to the main thread
    parentPort.postMessage(Papa.unparse([entry]) + '\n');
  }
} catch (err) {
  console.error('Worker error processing', path, err);
}
