import { rm } from 'fs/promises';

import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { DEFAULT_MAX_CONCURRENCY, FORCE, MAX_CONCURRENCY } from '../util/constants.js';
import { ensureDirsExist, fileExists } from '../util/fs.js';
import { logEvent } from '../util/logging.js';
import { DOWNLOAD_STEP, getDownloaderRef, getSummaryRef } from './utils.js';

/**
 * Attempts to download a dataset.
 * Saves the status of the download in the associated dataset summary.
 *
 * @param {Dataset} dataset Dataset to download
 */
async function tryDownload(dataset) {
  const force = dataset.config.get(FORCE, false);
  await ensureDirsExist(dataset.dirPath);
  // This precheck could be moved to before prepareDownloads to reduce excess work even further
  if (!force && (await fileExists(dataset.dataFilePath))) {
    await mergeExistingDataset(dataset);
    await dataset.save();
    markSuccess(dataset);
    return;
  }

  try {
    const downloader = getDownloaderRef(dataset);
    await logEvent('Download', dataset.id, () => downloader.download(dataset));
    if (!dataset.scratch.exclude) {
      await dataset.save();
    }
    markSuccess(dataset);
  } catch (error) {
    markFailure(dataset, error);
    await rm(dataset.dirPath, { recursive: true });
  }
}

/**
 * Marks the download step as successful
 *
 * @param {Dataset} dataset
 */
function markSuccess(dataset) {
  getSummaryRef(dataset).setSuccess(DOWNLOAD_STEP);
}

/**
 * Marks the download step as failed
 *
 * @param {Dataset} dataset
 * @param {any} error
 */
function markFailure(dataset, error) {
  getSummaryRef(dataset).setFailure(DOWNLOAD_STEP, error.message ?? error);
}

/**
 * Merges a dataset with values from file if the dataset file exists
 *
 * @param {Dataset} dataset
 */
async function mergeExistingDataset(dataset) {
  try {
    const existing = await Dataset.load(dataset.filePath, dataset.config);
    for (const key of Object.keys(existing)) {
      dataset[key] = dataset[key] || existing[key];
    }
  } catch {
    // Do nothing if there is no existing dataset.json file
  }
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
