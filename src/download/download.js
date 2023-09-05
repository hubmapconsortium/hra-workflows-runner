import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { fileExists } from '../util/fs.js';
import { DOWNLOAD_STEP, getDownloaderRef, getSummaryRef } from './utils.js';

async function tryDownload(dataset) {
  if (fileExists(dataset.dataFilePath)) {
    markSuccess(dataset);
    return;
  }

  try {
    const downloader = getDownloaderRef(dataset);
    await downloader.download(dataset);
    markSuccess(dataset);
  } catch (error) {
    markFailure(dataset, error);
  }
}

function markSuccess(dataset) {
  getSummaryRef(dataset).setSuccess(DOWNLOAD_STEP);
}

function markFailure(dataset, error) {
  const msg = `Download failed: ${error.message ?? error}`;
  getSummaryRef(dataset).setFailure(DOWNLOAD_STEP, msg);
}

/**
 * Downloads each dataset
 *
 * @param {Dataset} datasets
 * @param {Config} config
 */
export async function download(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  await concurrentMap(datasets, tryDownload, { maxConcurrency });
}
