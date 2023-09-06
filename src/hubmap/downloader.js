import { Config } from '../util/config.js';
import { FORCE } from '../util/constants.js';
import { downloadFile } from '../util/fs.js';
import { getMetadataLookup } from './metadata.js';

const HUBMAP_TOKEN = 'HUBMAP_TOKEN';
const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const HUBMAP_ASSETS_URL = 'HUBMAP_ASSETS_URL';
const DEFAULT_HUBMAP_SEARCH_URL =
  'https://search.api.hubmapconsortium.org/v3/portal/search';
const DEFAULT_HUBMAP_ASSETS_URL = 'https://assets.hubmapconsortium.org/';

export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.token = config.get(HUBMAP_TOKEN);
    /** @type {string} */
    this.searchUrl = config.get(HUBMAP_SEARCH_URL, DEFAULT_HUBMAP_SEARCH_URL);
    /** @type {string} */
    this.assetsUrl = config.get(HUBMAP_ASSETS_URL, DEFAULT_HUBMAP_ASSETS_URL);
  }

  async prepareDownload(datasets) {
    const ids = datasets.map((dataset) => dataset.id);
    const lookup = await getMetadataLookup(ids, this.searchUrl, this.token);
    for (const dataset of datasets) {
      const metadata = lookup.get(dataset.id);
      Object.assign(dataset, metadata);
    }
  }

  async download(dataset) {
    if (!dataset.uuid) {
      throw new Error('Missing uuid');
    }

    const url = new URL(`${dataset.uuid}/raw_expr.h5ad`, this.assetsUrl);
    url.searchParams.set('token', this.token);
    await downloadFile(dataset.dataFilePath, url, {
      overwrite: this.config.get(FORCE, false),
    });
  }
}
