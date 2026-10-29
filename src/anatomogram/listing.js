import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';
import { cacheExperiments, getExperimentFilePath } from './utils.js';

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.getListingScriptFile = 'get_dataset_listing.py';
    /** @type {string} */
    this.getListingScriptFilePath = getSrcFilePath(config, 'anatomogram', this.getListingScriptFile);
  }

  async getDatasets() {
    const config = this.config;
    await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));
    const experiments = await cacheExperiments(config);

    const datasets = await concurrentMap(experiments, async (code) => {
      const { stdout } = await execFile('python3', [
        this.getListingScriptFilePath,
        getExperimentFilePath(code, config),
        '--experiment',
        code,
      ]);
      return JSON.parse(stdout);
    });
    return datasets.flat().sort();
  }
}
