import { createScratchGetSet } from '../dataset/dataset.js';
import { DatasetSummary } from '../dataset/summary.js';
import { IDownloader } from '../util/handler.js';

export const DOWNLOAD_STEP = 'downloaded';
export const SUMMARY_REF = 'summary_ref';
export const DOWNLOADER_REF = 'downloader_ref';

export const { get: getSummaryRef, set: setSummaryRef } =
  /** @type {import('../dataset/dataset.js').ScratchGetSetPair<DatasetSummary>} */
  (createScratchGetSet(SUMMARY_REF));

export const { get: getDownloaderRef, set: setDownloaderRef } =
  /** @type {import('../dataset/dataset.js').ScratchGetSetPair<IDownloader>} */
  (createScratchGetSet(DOWNLOADER_REF));
