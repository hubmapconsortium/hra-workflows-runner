
export function supports(_dataset) {
  return true;
}

export class Downloader {
  async prepareDownload(_datasets) {
    return [];
  }

  async download(_dataset) {
    throw new Error('Not supported');
  }
}
