import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { getSrcFilePath } from '../util/paths.js';
import { Config } from '../util/config.js';
import { FORCE } from '../util/constants.js';
import { downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getMetadataLookup } from './metadata.js';

const HUBMAP_TOKEN = 'HUBMAP_TOKEN';
const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const HUBMAP_ASSETS_URL = 'HUBMAP_ASSETS_URL';
const DEFAULT_HUBMAP_SEARCH_URL =
  'https://search.api.hubmapconsortium.org/v3/portal/search';
const DEFAULT_HUBMAP_ASSETS_URL = 'https://assets.hubmapconsortium.org/';

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
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
    this.exprAdjustScript = 'expr_h5ad_adjust.py';
    /** @type {string} */
    this.exprAdjustScriptFilePath = getSrcFilePath(
      config,
      'hubmap',
      this.exprAdjustScript
    );
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
      throw new Error('Missing uuid - Dataset might have been deleted');
    }

    const url = new URL(`${dataset.uuid}/expr.h5ad`, this.assetsUrl);
    url.searchParams.set('token', this.token);
    await downloadFile(dataset.dataFilePath, url, {
      headers: {
        Authorization: `Bearer ${this.token}`,
      },
      overwrite: this.config.get(FORCE, false),
    });
    
    const { stdout } = await execFile('python3', [
      this.exprAdjustScriptFilePath,
      dataset.dataFilePath,
      '--assay',
      dataset.assay_type,
      '--output',
      dataset.dataFilePath,
    ]);
  }
  
}
