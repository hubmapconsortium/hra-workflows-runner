import { env } from 'node:process';
import Downloader from './downloader.js';
import { downloadRemoteFile } from './utils.js';

const HUBMAP_TOKEN = env['HUBMAP_TOKEN'];

export default class HubmapDownloader extends Downloader {
  static supports(dataset) {
    return /^hbm/i.test(dataset.id);
  }

  async doDownload() {
    const uuid = this.dataset.uuid;
    const url = `https://assets.hubmapconsortium.org/${uuid}/raw_expr.h5ad?token=${HUBMAP_TOKEN}`;
    if (!uuid) {
      throw new Error(`Could not find uuid for ${this.dataset.id}`);
    }

    await downloadRemoteFile(url, this.dataFile);
  }
}
