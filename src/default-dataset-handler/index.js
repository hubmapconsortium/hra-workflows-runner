import { IDownloader } from '../util/handler.js';

export function supports(_dataset) {
  return true;
}

/** @implements {IDownloader} */
export class Downloader {
  async prepareDownload(_datasets) {
    return [];
  }

  async download(_dataset) {
    throw new Error('Not supported');
  }
}
