import { execFile as callbackExecFile } from 'node:child_process';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { Config } from '../util/config.js';
import { FORCE } from '../util/constants.js';
import { downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';

const GTEX_FULL_DATA_URL = 'GTEX_FULL_DATA_URL';
const DEFAULT_GTEX_FULL_DATA_URL =
  'https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad';
const GTEX_DOI_URL = 'https://doi.org/10.1126/science.abl4290';

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class Downloader {
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
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(
      config,
      'gtex',
      this.extractScriptFile
    );
  }

  async prepareDownload(_datasets) {
    await downloadFile(this.dataFilePath, this.dataUrl, {
      overwrite: this.config.get(FORCE, false),
    });

    for (const dataset of _datasets) {
      const id = dataset.id;
      dataset.dataset_iri = `${GTEX_DOI_URL}#${id}`;
    }
  }

  async download(dataset) {
    const { stdout } = await execFile('python3', [
      this.extractScriptFilePath,
      this.dataFilePath,
      '--dataset',
      dataset.id,
      '--output',
      dataset.dataFilePath,
    ]);

    // Parse organ line. Format: `organ: X\n`
    const organ_match = /organ:(.+)\n/i.exec(stdout);
    dataset.organ = organ_match?.[1].trim() ?? '';

    // Parse sex line. Format: `sex: X\n`
    const sex_match = /sex:(.+)\n/i.exec(stdout);
    dataset.donor_sex = sex_match?.[1].trim() ?? '';

    // Parse age line. Format: `age: X\n`
    const age_match = /age:(.+)\n/i.exec(stdout);
    dataset.donor_age = age_match?.[1].trim() ?? '';
  }
}
