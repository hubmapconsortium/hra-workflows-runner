import { IDownloader, IJobGenerator, IListing } from '../util/handler.js';

export function supports(_dataset) {
  return true;
}

/** @implements {IListing} */
export class Listing {
  async getDatasets() {
    return [];
  }
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

/** @implements {IJobGenerator} */
export class JobGenerator {
  async prepareJobs(_datasets) {
    return [];
  }

  async createJob(_dataset) {
    throw new Error('Not supported');
  }
}
