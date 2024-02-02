import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { IDownloader } from '../util/handler.js';
import { diff, groupBy } from '../util/iter.js';
import { logEvent } from '../util/logging.js';
import { DOWNLOAD_STEP, getDownloaderRef, getSummaryRef } from './utils.js';

/**
 * Prepares datasets for downloading
 *
 * @param {[IDownloader, Dataset[]]} param0 Downloader and datasets packed into a tuple
 * @returns Datasets to download
 */
async function prepare([downloader, datasets]) {
  const newDatasets = await tryPrepare(downloader, datasets);
  if (newDatasets === undefined) {
    return datasets;
  } else if (newDatasets === 'error') {
    return [];
  }

  const notSupported = Array.from(diff(datasets, newDatasets));
  markNotSupported(notSupported);
  return newDatasets;
}

/**
 * Runs IDownloader#prepareDownload safely catching and handling any errors thrown
 *
 * @param {IDownloader} downloader Downloader instance
 * @param {Dataset[]} datasets Datasets
 * @returns Potentially filtered datasets or 'error' in case an error was caught
 */
async function tryPrepare(downloader, datasets) {
  try {
    return await logEvent(`PrepareDownload:${datasets[0].handler}`, () => downloader.prepareDownload(datasets));
  } catch (error) {
    markFailure(datasets, error);
    return 'error';
  }
}

/**
 * Marks multiple dataset's download step as not supported
 *
 * @param {Dataset[]} datasets Datasets
 */
function markNotSupported(datasets) {
  datasets.forEach((d) => getSummaryRef(d).setNotSupported(DOWNLOAD_STEP));
}

/**
 * Marks multiple dataset's download step as failed
 *
 * @param {Dataset[]} datasets Datasets
 * @param {any} error Cause of failure
 */
function markFailure(datasets, error) {
  const msg = `Prepare download failed: ${error.message ?? error}`;
  datasets.forEach((d) => getSummaryRef(d).setFailure(DOWNLOAD_STEP, msg));
}

/**
 * Prepares datasets for downloading
 *
 * @param {Dataset[]} datasets Datasets
 * @param {Config} config
 * @returns {Promise<Dataset[]>}
 */
export async function prepareDownloads(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const grouped = groupBy(datasets, getDownloaderRef);
  const downloadable = await concurrentMap(Array.from(grouped), prepare, {
    maxConcurrency,
  });
  return [].concat(...downloadable);
}
