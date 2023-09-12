import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from './concurrent-map.js';
import { Config } from './config.js';
import { DATASET_HANDLERS } from './constants.js';
import { getSrcFilePath } from './paths.js';

/**
 * Downloading interface implemented by handlers
 *
 * @interface
 */
export class IDownloader {
  /**
   * Prepares datasets for downloading
   *
   * @param {Dataset[]} datasets Datasets to prepare
   * @returns {Promise<Dataset[] | undefined>}
   */
  async prepareDownload(datasets) {}

  /**
   * Download a single dataset
   *
   * @param {Dataset} dataset Dataset to download
   * @throws If the download fails
   */
  async download(dataset) {}
}

/**
 * @typedef {Object} DatasetHandler Interface implemented by handlers
 * @property {function(Dataset): boolean} supports Test whether this handler supports a dataset
 * @property {function(new:IDownloader, Config)} Downloader Downloader class
 */

/**
 * Load dataset handler modules
 *
 * @param {Config} config Configuration
 * @returns {Promise<Map<string, DatasetHandler>>}
 */
export async function loadDatasetHandlers(config) {
  const handlers = await concurrentMap(DATASET_HANDLERS, (name) =>
    loadDatasetHandler(name, config)
  );
  return new Map(handlers.filter((h) => !!h));
}

async function loadDatasetHandler(name, config) {
  const path = getSrcFilePath(config, name, 'index.js');

  try {
    return [name, await import(path)];
  } catch {
    const msg = `Failed to load dataset handler '${name}' with path '${path}'`;
    console.warn(msg);
    return undefined;
  }
}
