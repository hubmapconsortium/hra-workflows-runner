import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { join } from 'node:path';
import { Config } from '../util/config.js';
import { ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { getDataRepoDir, getSrcFilePath } from '../util/paths.js';

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;

    /** @type {string} - Relative or absolute path to your map.json */
    this.mapFilePath = join(process.cwd(), 'disco_map.json');

    /** @type {string} */
    this.getListingScriptFile = 'get_dataset_listing.py';

    /** @type {string} */
    this.getListingScriptFilePath = getSrcFilePath(
      config,
      'disco',  
      this.getListingScriptFile
    );
  }

  async getDatasets() {
    const config = this.config;
    await ensureDirsExist(getDataRepoDir(config));  

    const { stdout } = await execFile('python3', [
      this.getListingScriptFilePath,
      this.mapFilePath,
    ]);

    return JSON.parse(stdout);
  }
}
