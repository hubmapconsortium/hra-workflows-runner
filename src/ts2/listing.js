import { execFile as callbackExecFile } from 'node:child_process';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';
import { cacheCollections } from './utils.js';

const TS2_FIGSHARE_ID = 'TS2_FIGSHARE_ID';
const DEFAULT_TS2_FIGSHARE_ID = '27921984';

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.figshareId = config.get(TS2_FIGSHARE_ID, DEFAULT_TS2_FIGSHARE_ID);
    /** @type {string} */
    this.getListingScriptFile = 'get_dataset_listing.py';
    /** @type {string} */
    this.getListingScriptFilePath = getSrcFilePath(config, 'ts2', this.getListingScriptFile);
  }

  async getDatasets() {
    const config = this.config;
    await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));
    const collections = await cacheCollections(this.figshareId, config);

    const datasets = await concurrentMap(collections, async (file) => {
      const dataFilePath = join(getCacheDir(config), 'ts2', file.name);
      const { stdout } = await execFile('python3', [this.getListingScriptFilePath, dataFilePath]);
      return JSON.parse(stdout);
    });
    return datasets.flat().sort();
  }
}
