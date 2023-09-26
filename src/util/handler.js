import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from './concurrent-map.js';
import { Config } from './config.js';
import { DATASET_HANDLERS, DEFAULT_DATASET_HANDLERS } from './constants.js';
import { getSrcFilePath } from './paths.js';

/**
 * @typedef {Object} DatasetHandler
 * Primary interface that must be exported by a handler's index.js entrypoint
 *
 * @property {function(Dataset): boolean} supports
 * Method to test whether this handler supports a specific dataset
 *
 * @property {function(new:IDownloader, Config)} Downloader
 * Class implementing the data downloading interface
 */

/**
 * Downloading interface implemented by handlers
 *
 * @interface
 */
export class IDownloader {
  /**
   * Prepares datasets for downloading.
   * This method can attach metadata about the dataset as well as
   * temporary data used to download the dataset's data to file.
   * It may also filter and return a subset of the passed
   * datasets should some of them not support downloading.
   * If the method throws all of the datasets will be marked as failed to download and
   * no calls to {@link download} will occur.
   *
   * @param {Dataset[]} datasets Datasets to prepare for downloading
   * @returns {Promise<Dataset[] | undefined>}
   * A subset of the initial datasets or undefined to download all datasets
   */
  async prepareDownload(datasets) {}

  /**
   * Download a single dataset to file.
   *
   * @param {Dataset} dataset Dataset to download
   * @throws If the download fails
   */
  async download(dataset) {}
}

/**
 * Load dataset handler modules
 *
 * @param {Config} config Configuration
 * @returns {Promise<Map<string, DatasetHandler>>}
 */
export async function loadDatasetHandlers(config) {
  const fallbackHandlerName = 'default-dataset-handler';
  const handlerNames = config.get(DATASET_HANDLERS, DEFAULT_DATASET_HANDLERS);
  const handlerNamesWithFallback = [...handlerNames, fallbackHandlerName];
  const handlerModules = await concurrentMap(handlerNamesWithFallback, (name) =>
    loadDatasetHandler(name, config)
  );
  const handlers = handlerModules.filter(verifyDatasetHandler);

  return new Map(handlers);
}

/**
 * Load a single handler module
 *
 * @param {string} name Handler name
 * @param {Config} config Configuration
 * @returns {Promise<[string, object] | undefined>}
 */
async function loadDatasetHandler(name, config) {
  const path = getSrcFilePath(config, name, 'index.js');

  try {
    return [name, await import(path)];
  } catch {
    const msg = `Failed to load dataset handler '${name}' using path '${path}'`;
    console.warn(msg);
    return undefined;
  }
}

/**
 * Tests whether a optional key/value pair is a dataset handler
 *
 * @param {[string, object] | undefined} maybeHandler Object to verify
 * @returns {maybeHandler is [string, DatasetHandler]}
 */
function verifyDatasetHandler(maybeHandler) {
  if (!maybeHandler) {
    return false;
  }

  const [name, module] = maybeHandler;
  const hasSupports = typeof module.supports === 'function';
  const hasDownloader = typeof module.Downloader === 'function';
  if (!hasSupports || !hasDownloader) {
    const msg = `Dataset handler '${name}' does not implement the required interfaces`;
    console.warn(msg);
    return false;
  }

  return true;
}
