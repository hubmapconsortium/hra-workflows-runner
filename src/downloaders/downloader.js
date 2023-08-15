import { mkdir, rm } from 'node:fs/promises';
import { join } from 'node:path';
import { fileExists } from './utils.js';

export default class Downloader {
  /**
   * Test whether this class supports downloading a dataset
   * @param {object} _dataset Dataset to test
   * @returns True if the dataset can be downloaded using this class, false otherwise
   */
  static supports(_dataset) {
    return false;
  }

  /**
   * Initialize the downloader
   * @param {object} dataset The dataset to download
   * @param {string} outDir Output directory
   * @param {string} cacheDir Temporary cache directory
   */
  constructor(dataset, outDir, cacheDir) {
    this.dataset = dataset;
    this.outDir = outDir;
    this.dataDir = join(outDir, dataset.id);
    this.dataFile = join(outDir, dataset.id, 'data.h5ad');
    this.cacheDir = cacheDir;
  }

  /**
   * Download the dataset passed in the constructor
   * @returns Resolves to undefined when the download is complete
   */
  async download() {
    if (await fileExists(this.dataFile)) {
      return;
    }

    try {
      await mkdir(this.dataDir, { recursive: true });
      await this.doDownload();
    } catch (error) {
      await rm(this.dataDir, { recursive: true, force: true });
      throw error;
    }
  }

  /**
   * Method to override in subclasses to perform the actual download
   */
  async doDownload() {
    throw new Error('doDownload() must be overridden in the subclass');
  }
}
