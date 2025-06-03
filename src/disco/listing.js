import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { join } from 'node:path';
import { Config } from '../util/config.js';
import { downloadFile, ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { FORCE } from '../util/constants.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';

const DISCO_METADATA_URL = "https://disco.bii.a-star.edu.sg/disco_v3_api/toolkit/getSampleMetadata";

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.cacheDir = getCacheDir(config);
    /** @type {string} */
    this.dataRepoDir = getDataRepoDir(config);
    /** @type {string} */
    this.metadataFile = "disco_metadata.tsv";
    /** @type {string} */
    this.metadataFilePath = join(this.cacheDir, this.metadataFile);
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
    // Ensure required directories exist
    await ensureDirsExist(this.dataRepoDir, this.cacheDir);

    // Download metadata file
    await downloadFile(this.metadataFilePath, DISCO_METADATA_URL, {
      overwrite: this.config.get(FORCE, false),
    });
    
    const args = [
      this.metadataFilePath,
      ...this.batchFiles
    ];
    const { stdout } = await execFile('python3', [
      this.getListingScriptFilePath,
      ...args
    ]);
    return JSON.parse(stdout);
  }
}
