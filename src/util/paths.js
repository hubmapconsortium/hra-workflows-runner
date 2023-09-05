import { join } from 'node:path';
import {
  ALGORITHM_REPORT_FILE,
  CACHE_DIR,
  DATA_FILE,
  DATA_REPO_DIR,
  DEFAULT_CACHE_DIR,
  LISTING_FILE,
  OUTPUT_DIR,
  SRC_DIR,
  SUMMARIES_FILE,
} from './constants.js';

export function getDirForId(id) {
  return id.replaceAll(/[^\w-._]/gi, '_');
}

export function getOutputDir(config) {
  return config.get(OUTPUT_DIR);
}

export function getDataRepoDir(config) {
  return config.get(DATA_REPO_DIR);
}

export function getDataDir(config, dir) {
  return join(getDataRepoDir(config), dir);
}

export function getCacheDir(config) {
  return config.get(CACHE_DIR, DEFAULT_CACHE_DIR);
}

export function getSrcDir(config) {
  return config.get(SRC_DIR);
}

export function getListingFilePath(config) {
  return join(getOutputDir(config), LISTING_FILE);
}

export function getSummariesFilePath(config) {
  return join(getOutputDir(config), SUMMARIES_FILE);
}

export function getDataFilePath(config, dir) {
  return join(getDataDir(config, dir), DATA_FILE);
}

export function getAlgorithmReportFilePath(config, dir, algorithm) {
  return join(getDataDir(config, dir), algorithm, ALGORITHM_REPORT_FILE);
}

export function getSrcFilePath(config, ...paths) {
  return join(getSrcDir(config), ...paths);
}
