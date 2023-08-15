import { execFile as rawExecFile } from 'node:child_process';
import { join } from 'node:path';
import { env } from 'node:process';
import { promisify } from 'node:util';
import Downloader from './downloader.js';
import { downloadRemoteFile, once } from './utils.js';

const FULL_DATA_URL =
  env['GTEX_FULL_DATA_URL'] ||
  'https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad';

const PYTHON_EXTRACT_SCRIPT = 'src/extract_gtex_dataset.py';

const execFile = promisify(rawExecFile);

const downloadPrimaryData = once((file) =>
  downloadRemoteFile(FULL_DATA_URL, file)
);

export default class GtexDownloader extends Downloader {
  static supports(dataset) {
    return /^gtex/i.test(dataset.id);
  }

  constructor(...args) {
    super(...args);
    this.fullDataFile = join(this.cacheDir, 'gtex-full-data.h5ad');
  }

  async doDownload() {
    await downloadPrimaryData(this.fullDataFile);
    await this.extractDataset();
  }

  async extractDataset() {
    return execFile('python3', [
      PYTHON_EXTRACT_SCRIPT,
      this.fullDataFile,
      '--dataset',
      this.dataset.id,
      '--output',
      this.dataFile,
    ]);
  }
}
