import { createScratchGetSet } from '../dataset/dataset.js';

export const DOWNLOAD_STEP = 'downloaded';
export const SUMMARY_REF = 'summary_ref';
export const DOWNLOADER_REF = 'downloader_ref';

export const { get: getSummaryRef, set: setSummaryRef } =
  createScratchGetSet(SUMMARY_REF);

export const { get: getDownloaderRef, set: setDownloaderRef } =
  createScratchGetSet(DOWNLOADER_REF);
