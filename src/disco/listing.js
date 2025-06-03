import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { Config } from '../util/config.js';
import { downloadFile, ensureDirsExist } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { FORCE } from '../util/constants.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';

const DISCO_BATCH_URLS = [
  "https://zenodo.org/records/14159931/files/batch_1.tar.gz?download=1",
  "https://zenodo.org/records/14160154/files/batch_2.tar.gz?download=1",
  "https://zenodo.org/records/14160213/files/batch_3.tar.gz?download=1",
  "https://zenodo.org/records/14160221/files/batch_4.tar.gz?download=1",
  "https://zenodo.org/records/14160748/files/batch_5.tar.gz?download=1",
  "https://zenodo.org/records/14160802/files/batch_6.tar.gz?download=1",
  "https://zenodo.org/records/14166702/files/batch_7.tar.gz?download=1",
  "https://zenodo.org/records/15236185/files/batch_8.tar.gz?download=1",
  "https://zenodo.org/records/15236615/files/batch_9.tar.gz?download=1"
];

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
    /** @type {string[]} */
    this.batchFiles = DISCO_BATCH_URLS.map(url => {
      const batchName = url.match(/batch_\d+\.tar\.gz/)[0];
      return join(this.cacheDir, batchName);
    });
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

    // Download all batch tar files
    for (let i = 0; i < DISCO_BATCH_URLS.length; i++) {
      await downloadFile(this.batchFiles[i], DISCO_BATCH_URLS[i], {
        overwrite: this.config.get(FORCE, false),
      });
    }

    // Call the Python script to extract the listing
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
