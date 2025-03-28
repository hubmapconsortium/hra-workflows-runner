import { join } from 'node:path';
import { Config } from './config.js';
import {
  ALGORITHM_ANNOTATIONS_FILE,
  ALGORITHM_REPORT_FILE,
  ALGORITHM_SUMMARY_JSON_LD_FILE,
  CACHE_DIR,
  CROSSWALKING_TABLES_DIR,
  DATASET,
  DATASETS_DIR,
  DATASET_FILE,
  DATASET_LIST,
  DATA_FILE,
  DATA_REPO_DIR,
  DEFAULT_CACHE_DIR,
  DEFAULT_DATASET_LIST,
  JOB_FILE,
  LISTING_FILE,
  MODELS_DIR,
  OUTPUT_DIR,
  SRC_DIR,
  SUMMARIES_FILE,
} from './constants.js';

/**
 * Turns a dataset identifier into a directory
 *
 * @param {string} id Dataset identifier
 */
export function getDirForId(id) {
  return id.replaceAll(/[^\w-._]/gi, '_');
}

/**
 * Get the configured dataset directory
 *
 * @param {Config} config Configuration
 */
export function getDatasetDir(config) {
  return join(config.get(DATASETS_DIR), config.get(DATASET));
}

/**
 * Get the configured output directory
 *
 * @param {Config} config Configuration
 */
export function getOutputDir(config) {
  return config.get(OUTPUT_DIR);
}

/**
 * Get the configured data repository directory
 *
 * @param {Config} config Configuration
 */
export function getDataRepoDir(config) {
  return config.get(DATA_REPO_DIR);
}

/**
 * Get the configured data directory
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 */
export function getDataDir(config, dir) {
  return join(getDataRepoDir(config), dir);
}

/**
 * Get the configured cache directory
 *
 * @param {Config} config Configuration
 */
export function getCacheDir(config) {
  return config.get(CACHE_DIR, DEFAULT_CACHE_DIR);
}

/**
 * Get the configured models directory
 *
 * @param {Config} config Configuration
 */
export function getModelsDir(config) {
  return config.get(MODELS_DIR);
}

/**
 * Get the configured source directory
 *
 * @param {Config} config Configuration
 */
export function getSrcDir(config) {
  return config.get(SRC_DIR);
}

/**
 * Get the configured list file path
 *
 * @param {Config} config Configuration
 */
export function getDatasetListFilePath(config) {
  return join(getDatasetDir(config), config.get(DATASET_LIST, DEFAULT_DATASET_LIST));
}

/**
 * Get the configured listing file path
 *
 * @param {Config} config Configuration
 */
export function getListingFilePath(config) {
  return join(getOutputDir(config), LISTING_FILE);
}

/**
 * Get the configured summaries file path
 *
 * @param {Config} config Configuration
 */
export function getSummariesFilePath(config) {
  return join(getOutputDir(config), SUMMARIES_FILE);
}

/**
 * Get the configured dataset file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 */
export function getDatasetFilePath(config, dir) {
  return join(getDataDir(config, dir), DATASET_FILE);
}

/**
 * Get the configured data file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 */
export function getDataFilePath(config, dir) {
  return join(getDataDir(config, dir), DATA_FILE);
}

/**
 * Get the configured job file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 */
export function getJobFilePath(config, dir) {
  return join(getDataDir(config, dir), JOB_FILE);
}

/**
 * Get the configured job file path with added suffix
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 * @param {string} suffix Suffix
 */
export function getJobFilePathWithSuffix(config, dir, suffix) {
  const [name, ext] = JOB_FILE.split('.', 2);
  const file = `${name}-${suffix}.${ext}`;
  return join(getDataDir(config, dir), file);
}

/**
 * Get the configured algorithm report file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 * @param {string} algorithm Algorithm
 */
export function getAlgorithmReportFilePath(config, dir, algorithm) {
  return join(getDataDir(config, dir), algorithm, ALGORITHM_REPORT_FILE);
}

/**
 * Get the configured algorithm summary file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 * @param {string} algorithm Algorithm
 */
export function getAlgorithmSummaryJsonLdFilePath(config, dir, algorithm) {
  return join(getDataDir(config, dir), algorithm, ALGORITHM_SUMMARY_JSON_LD_FILE);
}

/**
 * Get the configured algorithm annotations file path
 *
 * @param {Config} config Configuration
 * @param {string} dir Dataset directory
 * @param {string} algorithm Algorithm
 */
export function getAlgorithmAnnotationsFilePath(config, dir, algorithm) {
  return join(getDataDir(config, dir), algorithm, ALGORITHM_ANNOTATIONS_FILE);
}

/**
 * Get the path to a source file
 *
 * @param {Config} config Configuration
 * @param {...string} paths Additional path segments
 */
export function getSrcFilePath(config, ...paths) {
  return join(getSrcDir(config), ...paths);
}

/**
 * Get the Crosswalking Table file
 *
 * @param {Config} config Configuration
 * @param {string} algorithm Algorithm
 */
export function getCrosswalkingFilePath(config, algorithm) {
  return join(config.get(CROSSWALKING_TABLES_DIR), `${algorithm}.csv`);
}
