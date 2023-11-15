import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { Config } from '../util/config.js';
import { downloadFile } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { FORCE } from '../util/constants.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';

const GTEX_FULL_DATA_URL = 'GTEX_FULL_DATA_URL';
const DEFAULT_GTEX_FULL_DATA_URL =
  'https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad';

const execFile = promisify(callbackExecFile);

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.dataUrl = config.get(GTEX_FULL_DATA_URL, DEFAULT_GTEX_FULL_DATA_URL);
    /** @type {string} */
    this.dataUrlHash = createHash('md5').update(this.dataUrl).digest('hex');
    /** @type {string} */
    this.dataFile = `gtex-full-data-${this.dataUrlHash}.h5ad`;
    /** @type {string} */
    this.dataFilePath = join(getCacheDir(config), this.dataFile);
    /** @type {string} */
    this.getListingScriptFile = 'get_dataset_listing.py';
    /** @type {string} */
    this.getListingScriptFilePath = getSrcFilePath(
      config,
      'gtex',
      this.getListingScriptFile
    );
  }

  async getDatasets() {
    await downloadFile(this.dataFilePath, this.dataUrl, {
      overwrite: this.config.get(FORCE, false),
    });
    const { stdout } = await execFile('python3', [
      this.getListingScriptFilePath,
      this.dataFilePath,
    ]);
    return JSON.parse(stdout);
  }
}
