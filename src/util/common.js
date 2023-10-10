import { readFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { Dataset, createScratchGetSet } from '../dataset/dataset.js';
import { DatasetSummaries, DatasetSummary } from '../dataset/summary.js';
import { concurrentMap } from './concurrent-map.js';
import { Config } from './config.js';
import {
  DEFAULT_MAX_CONCURRENCY,
  MAX_CONCURRENCY,
  REQUIRED_ENV_VARIABLES,
} from './constants.js';
import { defaultEnvReviver } from './default-env-reviver.js';
import {
  getDatasetFilePath,
  getDirForId,
  getListingFilePath,
  getSummariesFilePath,
} from './paths.js';

/**
 * Gets a new configuration object with the environment loaded and validated
 */
export function getConfig() {
  return new Config()
    .loadEnv(defaultEnvReviver)
    .validate(REQUIRED_ENV_VARIABLES);
}

export const { get: getSummaryRef, set: setSummaryRef } =
  /** @type {import('../dataset/dataset.js').ScratchGetSetPair<DatasetSummary | undefined>} */
  (createScratchGetSet('summaryRef'));

/**
 * Load a csv file
 *
 * @template T
 * @param {import('node:fs').PathLike} path Path to csv file
 * @returns {Promise<T[]>}
 */
export async function loadCsv(path) {
  const content = await readFile(path, { encoding: 'utf8' });
  const { data } = Papa.parse(content, { header: true, skipEmptyLines: 'greedy' });
  return data;
}

/**
 * Load a json file
 *
 * @template T
 * @param {import('node:fs').PathLike} path Path to json file
 * @returns {Promise<T>}
 */
export async function loadJson(path) {
  const content = await readFile(path, { encoding: 'utf8' });
  return JSON.parse(content);
}

/**
 * Load the dataset's listing file
 *
 * @param {Config} config Configuration
 * @returns {Promise<{ id: string; [key: string]: string }[]>}
 */
export async function loadListing(config) {
  const prefix = 'DATASET_COLUMN_';
  const listing = await loadCsv(getListingFilePath(config));
  const keys = config.getPrefixedKeys(prefix);
  const columns = keys.map((key) => config.get(key));
  const props = keys.map((key) => key.slice(prefix.length).toLowerCase());

  return listing.map((row) => parseListingRow(row, columns, props));
}

function parseListingRow(row, columns, props) {
  return columns.reduce(
    (res, col, index) => ({
      ...res,
      [props[index]]: row[col],
    }),
    {}
  );
}

/**
 * Load summaries from file
 *
 * @param {Config} config Configuration
 */
export async function loadSummaries(config) {
  const path = getSummariesFilePath(config);
  return DatasetSummaries.load(path);
}

/**
 * Save summaries to file
 *
 * @param {DatasetSummaries} summaries New summaries
 * @param {Config} config Configuration
 */
export async function saveSummaries(summaries, config) {
  const path = getSummariesFilePath(config);
  await summaries.save(path);
}

/**
 * Loads multiple datasets from file for the specified directories
 * or by id for each summary
 *
 * @param {string[] | DatasetSummary[] | DatasetSummaries} directoriesOrSummaries
 * @param {Config} config Configuration
 * @returns {Promise<Dataset[]>}
 */
export async function loadDatasets(directoriesOrSummaries, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const items = Array.isArray(directoriesOrSummaries)
    ? directoriesOrSummaries
    : Array.from(directoriesOrSummaries.values());

  return concurrentMap(items, (item) => loadDataset(item, config), {
    maxConcurrency,
  });
}

/**
 * Saves multiple datasets to file
 *
 * @param {Dataset[]} datasets Datasets
 * @param {Config} config Configuration
 */
export async function saveDatasets(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  await concurrentMap(datasets, (dataset) => dataset.save(), {
    maxConcurrency,
  });
}

async function loadDataset(dirOrSummary, config) {
  const [summary, directory] =
    typeof dirOrSummary === 'string'
      ? [undefined, dirOrSummary]
      : [dirOrSummary, getDirForId(dirOrSummary.id)];
  const path = getDatasetFilePath(config, directory);
  const dataset = await Dataset.load(path, config);

  setSummaryRef(dataset, summary);
  return dataset;
}
