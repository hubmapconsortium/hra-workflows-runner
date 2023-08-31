import { readFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { Config } from './config.js';
import {
  DATASET_ID_COLUMN,
  DATASET_LINK_COLUMN,
  REQUIRED_ENV_VARIABLES,
} from './constants.js';
import { defaultEnvReviver } from './env.js';
import { getListingFilePath } from './paths.js';

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
  const listing = await loadCsv(getListingFilePath(config));
  const idColumn = config.get(DATASET_ID_COLUMN);
  const linkColumn = config.get(DATASET_LINK_COLUMN);
  return listing.map((row) => ({ id: row[idColumn], link: row[linkColumn] }));
}
