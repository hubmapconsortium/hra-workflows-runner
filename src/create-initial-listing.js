import { writeFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { getConfig } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { Config } from './util/config.js';
import { DATASET_COLUMN_ID, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from './util/constants.js';
import { ensureDirsExist } from './util/fs.js';
import { loadDatasetHandlers } from './util/handler.js';
import { getDatasetListFilePath, getOutputDir } from './util/paths.js';

/**
 * Get all datasets listed by a handler
 *
 * @param {[string, import('./util/handler.js').DatasetHandler]} param0
 * @param {Config} config
 */
async function getDatasets([name, handler], config) {
  try {
    const listing = new handler.Listing(config);
    return await listing.getDatasets();
  } catch (error) {
    console.error(`Failed to get datasets with handler ${name}`, error);
    return [];
  }
}

async function main() {
  const config = getConfig().validate([DATASET_COLUMN_ID]);
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  await ensureDirsExist(getOutputDir(config));

  const handlers = await loadDatasetHandlers(config);
  const handlerEntries = Array.from(handlers.entries());
  const datasets = await concurrentMap(handlerEntries, (entry) => getDatasets(entry, config), {
    maxConcurrency,
  });

  const columnName = config.get(DATASET_COLUMN_ID);
  const ids = datasets.flat().map((id) => ({ [columnName]: id }));
  const data = Papa.unparse({
    data: ids,
    fields: [columnName],
  });

  const listFilePath = getDatasetListFilePath(config);
  await writeFile(listFilePath, data);
}

main();
