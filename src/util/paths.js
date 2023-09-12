import { join } from 'node:path';
import { Config } from './config.js';
import {
  ALGORITHM_REPORT_FILE,
  CACHE_DIR,
  DATASET_FILE,
  DATA_FILE,
  DATA_REPO_DIR,
  DEFAULT_CACHE_DIR,
  LISTING_FILE,
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
 * Get the configured source directory
 *
 * @param {Config} config Configuration
 */
export function getSrcDir(config) {
  return config.get(SRC_DIR);
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
 * Get the path to a source file
 *
 * @param {Config} config Configuration
 * @param {...string} paths Additional path segments
 */
export function getSrcFilePath(config, ...paths) {
  return join(getSrcDir(config), ...paths);
}
