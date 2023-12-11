import { Dataset } from '../dataset/dataset.js';
import { DatasetSummary } from '../dataset/summary.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { IJobGenerator } from '../util/handler.js';
import { groupBy } from '../util/iter.js';
import { getJobGeneratorRef, getSummaryRefMany } from './utils.js';

/**
 * Prepares datasets for creating job file definitions
 *
 * @param {[IJobGenerator, Dataset[]]} param0 Job generator and datasets packed into a tuple
 * @returns Prepared datasets
 */
async function prepare([generator, datasets]) {
  const newDatasets = await tryPrepare(generator, datasets);
  if (newDatasets === undefined) {
    return datasets;
  }

  const notSupported = Array.from(diff(datasets, newDatasets));
  const refs = getSummaryRefMany(notSupported);
  ALGORITHMS.forEach((step) => DatasetSummary.setNotSupportedMany(refs, step));

  return newDatasets;
}

/**
 * Prepare datasets while catching and handling errors
 *
 * @param {IJobGenerator} generator Job generator instance
 * @param {Dataset[]} datasets Datasets
 * @returns Prepared datasets or undefined on error
 */
async function tryPrepare(generator, datasets) {
  try {
    return await generator.prepareJobs(datasets);
  } catch (error) {
    const refs = getSummaryRefMany(datasets);
    ALGORITHMS.forEach((step) => DatasetSummary.setFailureMany(refs, step, error.message ?? error));
    return undefined;
  }
}

/**
 * Prepares datasets for job generation
 *
 * @param {Dataset[]} datasets
 * @param {Config} config
 * @returns {Promise<Dataset[]>}
 */
export async function prepareJobs(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const grouped = groupBy(datasets, getJobGeneratorRef);
  const runnable = await concurrentMap(Array.from(grouped), prepare, {
    maxConcurrency,
  });
  return [].concat(...runnable);
}
