import { rm } from 'fs/promises';
import { Dataset } from '../dataset/dataset.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import {
  DEFAULT_MAX_CONCURRENCY,
  FORCE,
  MAX_CONCURRENCY,
} from '../util/constants.js';
import { ensureDirsExist, fileExists } from '../util/fs.js';
import { logEvent } from '../util/logging.js';
import { DOWNLOAD_STEP, getDownloaderRef, getSummaryRef } from './utils.js';

async function tryDownload(dataset) {
  const force = dataset.config.get(FORCE, false);
  await ensureDirsExist(dataset.dirPath);
  // This precheck could be moved to before prepareDownloads to reduce excess work even further
  if (!force && (await fileExists(dataset.dataFilePath))) {
    if (!(await fileExists(dataset.filePath))) {
      await dataset.save();
    }

    markSuccess(dataset);
    return;
  }

  try {
    const downloader = getDownloaderRef(dataset);
    await logEvent('Download', dataset.id, () => downloader.download(dataset));
    await dataset.save();
    markSuccess(dataset);
  } catch (error) {
    markFailure(dataset, error);
    await rm(dataset.dirPath, { recursive: true });
  }
}

function markSuccess(dataset) {
  getSummaryRef(dataset).setSuccess(DOWNLOAD_STEP);
}

function markFailure(dataset, error) {
  getSummaryRef(dataset).setFailure(DOWNLOAD_STEP, error.message ?? error);
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
