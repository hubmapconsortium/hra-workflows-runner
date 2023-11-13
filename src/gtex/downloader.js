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

const ORGAN_MAPPING = {
  bladder: 'UBERON:0001255',
  blood: 'UBERON:0000178',
  bone_marrow: 'UBERON:0002371',
  eye: 'UBERON:0000970',
  heart: 'UBERON:0000948',
  large_intestine: 'UBERON:0000059',
  liver: 'UBERON:0002107',
  lung: 'UBERON:0002048',
  lymph_node: 'UBERON:0000029',
  mammary: 'UBERON:0001911',
  pancreas: 'UBERON:0001264',
  prostate: 'UBERON:0002367',
  skin: 'UBERON:0002097',
  small_intestine: 'UBERON:0002108',
  spleen: 'UBERON:0002106',
  thymus: 'UBERON:0002370',
  trachea: 'UBERON:0003126',
  uterus: 'UBERON:0000995',
  vasculature: 'UBERON:0004537',
};

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
    dataset.organ_source = organ_match?.[1].trim() ?? '';

    const organ = dataset.organ_source.toLowerCase();
    const uberon = ORGAN_MAPPING[organ];
    if (uberon === undefined) {
      const msg = `Unknown organ code '${dataset.organ}'`;
      throw new Error(msg);
    }
    dataset.organ = uberon;

    // Parse sex line. Format: `sex: X\n`
    const sex_match = /sex:(.+)\n/i.exec(stdout);
    dataset.donor_sex = sex_match?.[1].trim() ?? '';

    // Parse age line. Format: `age: X\n`
    const age_match = /age:(.+)\n/i.exec(stdout);
    dataset.donor_age = age_match?.[1].trim() ?? '';
  }
}
