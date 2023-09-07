import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { diff, groupBy } from '../util/iter.js';
import { DOWNLOAD_STEP, getDownloaderRef, getSummaryRef } from './utils.js';

async function prepare([downloader, datasets]) {
  const newDatasets = await tryPrepare(downloader, datasets);
  if (newDatasets === undefined) {
    return datasets;
  }

  const notSupported = Array.from(diff(datasets, newDatasets));
  markNotSupported(notSupported);
  return newDatasets;
}

async function tryPrepare(downloader, datasets) {
  try {
    return await downloader.prepareDownload(datasets);
  } catch (error) {
    markFailure(datasets, error);
    return undefined;
  }
}

function markNotSupported(datasets) {
  datasets.forEach((d) => getSummaryRef(d).setNotSupported(DOWNLOAD_STEP));
}

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
