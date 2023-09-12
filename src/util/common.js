import { readFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { Config } from './config.js';
import { REQUIRED_ENV_VARIABLES } from './constants.js';
import { defaultEnvReviver } from './default-env-reviver.js';
import { getListingFilePath } from './paths.js';

/**
 * Gets a new configuration object with the environment loaded and validated
 */
export function getConfig() {
  return new Config()
    .loadEnv(defaultEnvReviver)
    .validate(REQUIRED_ENV_VARIABLES);
}

/**
 * Load a csv file
 *
 * @template T
 * @param {import('node:fs').PathLike} path Path to csv file
 * @returns {Promise<T[]>}
 */
export async function loadCsv(path) {
  const content = await readFile(path, { encoding: 'utf8' });
  const { data } = Papa.parse(content, { header: true });
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
