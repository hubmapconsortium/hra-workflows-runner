import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { Config } from '../util/config.js';
import { FORCE } from '../util/constants.js';
import { downloadFile } from '../util/fs.js';
import { getSrcFilePath } from '../util/paths.js';

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class XConsortiaDownloader {
  constructor(config, token, searchUrl, assetsUrl) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.token = token;
    /** @type {string} */
    this.searchUrl = searchUrl;
    /** @type {string} */
    this.assetsUrl = assetsUrl;
    /** @type {string} */
    this.exprAdjustScript = 'expr_h5ad_adjust.py';
    /** @type {string} */
    this.exprAdjustScriptFilePath = getSrcFilePath(config, 'xconsortia', this.exprAdjustScript);
  }

  async prepareDownload(datasets) {
    const ids = datasets.map((dataset) => dataset.id);
    const lookup = await this.getMetadataLookup(ids);
    for (const dataset of datasets) {
      const metadata = lookup.get(dataset.id);
      Object.assign(dataset, metadata);
    }
    return datasets;
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

    await execFile('python3', [
      this.exprAdjustScriptFilePath,
      dataset.dataFilePath,
      '--assay',
      dataset.assay_type,
      '--output',
      dataset.dataFilePath,
    ]);
  }

  async getMetadataLookup(_ids) {
    throw new Error('getMetadataLookup must be overriden');
  }
}
