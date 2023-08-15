import Downloader from './downloader.js';

export default class CellXGeneDownloader extends Downloader {
  static supports(dataset) {
    return false;
  }

  async doDownload() {
    //
  }
}