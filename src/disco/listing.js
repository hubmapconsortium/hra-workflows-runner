import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { Config } from '../util/config.js';
import { ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { FORCE } from '../util/constants.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';

const DISCO_DATA_MAPPING = 'DISCO_DATA_MAPPING';
const DEFAULT_DISCO_DATA_MAPPING = 'path/to/mapping.json'; // This should be updated with actual mapping file path

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.mappingUrl = config.get(DISCO_DATA_MAPPING, DEFAULT_DISCO_DATA_MAPPING);
    /** @type {string} */
    this.mappingUrlHash = createHash('md5').update(this.mappingUrl).digest('hex');
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
    await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));
    
    // Execute Python script that reads the mapping file and returns dataset information
    const { stdout } = await execFile('python3', [
      this.getListingScriptFilePath,
      this.mappingUrl,
    ]);
    return JSON.parse(stdout);
  }
}
