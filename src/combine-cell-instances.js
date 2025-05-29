import { createWriteStream, existsSync } from 'node:fs';
import { join } from 'node:path';
import { createGzip } from 'node:zlib';
import Papa from 'papaparse';
import { DatasetSummaries } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { Config } from './util/config.js';
import { ALGORITHMS, POPV_METHOD_COUNT } from './util/constants.js';
import { readCsv } from './util/csv.js';
import {
  getAlgorithmAnnotationsFilePath,
  getDatasetFilePath,
  getDirForId,
  getOutputDir,
  getSummariesFilePath,
} from './util/paths.js';

/**
 * Attempts to read the summary jsonld file for a specific dataset and algorithm
 *
 * @param {string} dir Dataset directory
 * @param {string} algorithm Algorithm name
 * @param {Config} config Configuration
 * @returns The an async iterator summary file content or undefined if the file does not exist
 */
function readCellAnnotations(dir, algorithm, config) {
  const path = getAlgorithmAnnotationsFilePath(config, dir, algorithm);
  try {
    if (existsSync(path)) {
      return readCsv(path);
    }
    return undefined;
  } catch (err) {
    console.log('Error processing', path, err);
    return undefined;
  }
}

/**
 * Attempts to get metadata for the given dataset
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 * @returns the dataset IRI
 */
async function getDataset(dir, config) {
  const path = getDatasetFilePath(config, dir);
  try {
    return await loadJson(path);
  } catch {
    return undefined;
  }
}

async function main() {
  const config = getConfig();
  const popvMethodCount = parseInt(config.get(POPV_METHOD_COUNT, '6'));
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const outputFile = join(getOutputDir(config), 'sc-transcriptomics-cell-instances.csv.gz');
  let output = createWriteStream(outputFile, { autoClose: true });
  if (outputFile.endsWith('.gz')) {
    const gzip = createGzip();
    gzip.pipe(output);
    output = gzip;
  }
  const header = ['dataset', 'organ', 'tool', 'cell', 'cell_id', 'cell_label', 'match_type', 'confidence_score'];
  output.write(Papa.unparse([header]) + '\n');
  for (const item of summaries.values()) {
    const directory = getDirForId(item.id);
    const dataset = await getDataset(directory, config);
    for (const tool of ALGORITHMS) {
      const cells = readCellAnnotations(directory, tool, config) ?? [];
      for await (const cell of cells) {
        const confidence_score =
          parseFloat(cell['mapping.score']) ||
          parseFloat(cell.conf_score) ||
          parseFloat(cell.popv_prediction_score) / popvMethodCount ||
          0;
        const entry = {
          dataset: dataset.dataset_id,
          organ: dataset.organ ?? '',
          tool,
          cell: cell[''],
          cell_id: cell.clid,
          cell_label: cell.hra_prediction,
          match_type: cell.match_type,
          confidence_score,
        };
        const content = Papa.unparse([header.map((field) => entry[field])]);
        output.write(content + '\n');
      }
      if (output.flush) {
        // Drain gzip buffer periodically
        await new Promise((resolve) => {
          output.flush(resolve);
        });
      } else {
        // Drain buffer periodically
        await new Promise((resolve) => {
          output.once('drain', resolve);
        });
      }
    }
  }
  output.end();
}

main();
