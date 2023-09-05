import { readFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { concurrentMap } from './concurrent-map.js';
import { Config } from './config.js';
import { DATASET_HANDLERS, REQUIRED_ENV_VARIABLES } from './constants.js';
import * as DefaultDatasetHandler from './default-dataset-handler.js';
import { defaultEnvReviver } from './default-env-reviver.js';
import { getListingFilePath, getSrcFilePath } from './paths.js';

export function getConfig() {
  return new Config()
    .loadEnv(defaultEnvReviver)
    .validate(REQUIRED_ENV_VARIABLES);
}

export async function loadCsv(path) {
  const content = await readFile(path, { encoding: 'utf8' });
  const { data } = Papa.parse(content, { header: true });
  return data;
}

export async function loadJson(path) {
  const content = await readFile(path, { encoding: 'utf8' });
  return JSON.parse(content);
}

export async function loadListing(config) {
  const prefix = 'DATASET_COLUMN_';
  const listing = await loadCsv(getListingFilePath(config));
  const keys = config.getPrefixedKeys(prefix);
  const columns = keys.map((key) => config.get(key));
  const props = keys.map((key) => key.slice(prefix.length).toLowerCase());

  return listing.map((row) => parseListingRow(row, columns, props));
}

function parseListingRow(row, columns, props) {
  return columns.reduce((res, col, index) => ({
    ...res,
    [props[index]]: row[col],
  }), {});
}

export const DEFAULT_DATASET_HANDLER_NAME = '<default-dataset-handler>';

/**
 * Load dataset handler modules
 *
 * @param {Config} config 
 * @returns {Promise<Map<string, any>>}
 */
export async function loadDatasetHandlers(config) {
  const handlers = await concurrentMap(DATASET_HANDLERS, (name) =>
    loadDatasetHandler(name, config)
  );
  return new Map([
    ...handlers.filter((h) => !!h),
    [DEFAULT_DATASET_HANDLER_NAME, DefaultDatasetHandler],
  ]);
}

async function loadDatasetHandler(name, config) {
  try {
    const path = getSrcFilePath(config, name, 'index.js');
    return [name, await import(path)];
  } catch {
    return undefined;
  }
}
