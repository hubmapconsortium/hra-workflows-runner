/**
 * Entry point for combining single-cell transcriptomics annotations into a single CSV.
 *
 * This script:
 * - Loads dataset summaries
 * - Launches worker threads (via Node.js Worker API) to read and process per-tool annotations
 * - Writes out the results to a compressed or plain CSV stream
 * - Controls concurrency using a configurable MAX_PROCESSES limit
 * - Tracks and prints progress to stdout
 *
 * Environment Variables:
 * - MAX_PROCESSES: Maximum number of worker threads to run concurrently (default: all available cores)
 */

import { createWriteStream, existsSync } from 'node:fs';
import { cpus } from 'node:os';
import { join } from 'node:path';
import { Worker } from 'node:worker_threads';
import { createGzip } from 'node:zlib';
import Papa from 'papaparse';
import { DatasetSummaries } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { ALGORITHMS, POPV_METHOD_COUNT } from './util/constants.js';
import {
  getAlgorithmAnnotationsFilePath,
  getDatasetFilePath,
  getDirForId,
  getOutputDir,
  getSummariesFilePath,
} from './util/paths.js';

/**
 * Attempts to load dataset metadata as a JSON object.
 * @param {string} dir - Directory of the dataset
 * @param {object} config - Config object
 * @returns {Promise<object|undefined>}
 */
async function getDataset(dir, config) {
  try {
    return await loadJson(getDatasetFilePath(config, dir));
  } catch {
    return undefined;
  }
}

/**
 * Runs a single worker thread for a given task, wiring up event listeners to receive CSV lines.
 * @param {object} task - Task data including dataset info and onLine callback
 * @returns {Promise<void>}
 */
function runWorker(task, onMessage) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./combine-cell-instances.worker.js', import.meta.url), {
      workerData: task,
    });
    worker.on('message', onMessage);
    worker.on('error', reject);
    worker.on('exit', (code) => {
      if (code !== 0) {
        reject(new Error(`Worker exited with code ${code}`));
      } else {
        resolve();
      }
      worker.terminate();
    });
  });
}

/**
 * Executes a pool of async tasks with controlled concurrency.
 * @param {object[]} tasks - List of task objects to run
 * @param {number} maxConcurrent - Maximum number of concurrent workers
 * @returns {Promise<void>}
 */
async function runWithPool(tasks, onMessage, maxConcurrent = 10) {
  const queue = [...tasks];
  let completed = 0;

  // Prints dynamic progress to the same terminal line
  function printProgress() {
    process.stdout.write(`\rProcessed ${completed} / ${tasks.length} tasks`);
  }

  /**
   * Starts the next task if available, maintaining the pool limit.
   */
  async function taskRunner() {
    while (queue.length > 0) {
      const task = queue.shift();
      try {
        await runWorker(task, onMessage);
      } catch (err) {
        console.log(err);
      } finally {
        completed++;
        printProgress();
      }
    }
  }

  await Promise.all(Array.from({ length: Math.min(maxConcurrent, queue.length) }, taskRunner));
  process.stdout.write('\n');
}

/**
 * Main script logic: loads summaries, builds tasks, and coordinates worker processing.
 */
async function main() {
  const config = getConfig();
  const maxProcesses = parseInt(process.env.MAX_PROCESSES || cpus().length.toString(), 10);
  const popvMethodCount = parseInt(config.get(POPV_METHOD_COUNT, '6'));

  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const outputFile = join(getOutputDir(config), 'sc-transcriptomics-cell-instances2.csv.gz');

  // Create output stream (optionally gzipped)
  let output = createWriteStream(outputFile);
  if (outputFile.endsWith('.gz')) {
    const gzip = createGzip();
    gzip.pipe(output);
    output = gzip;
  }
  output.setMaxListeners(maxProcesses * 3);

  const tasks = [];

  // Generate all (dataset, tool) processing jobs
  for (const item of summaries.values()) {
    const directory = getDirForId(item.id);
    const dataset = await getDataset(directory, config);
    if (!dataset) continue;

    for (const tool of ALGORITHMS) {
      const path = getAlgorithmAnnotationsFilePath(config, directory, tool);
      if (existsSync(path)) {
        tasks.push({
          directory,
          tool,
          dataset,
          popvMethodCount,
          path,
        });
      }
    }
  }

  const header = ['dataset', 'organ', 'tool', 'cell', 'cell_id', 'cell_label', 'match_type', 'confidence_score'];
  output.write(Papa.unparse([header]) + '\n');

  // Poor man's shuffle
  tasks.sort(() => Math.random() - 0.5);

  console.log(`Starting ${tasks.length} tasks with max concurrency ${maxProcesses}...`);
  await runWithPool(tasks, (line) => output.write(line), maxProcesses);
  output.end();
  console.log('Done!');
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
