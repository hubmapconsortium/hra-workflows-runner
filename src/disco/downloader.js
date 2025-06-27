import { execFile as callbackExecFile, exec as rawExec } from 'node:child_process';
import fs, { existsSync } from 'node:fs';
import { join } from 'node:path';
import { promisify } from 'node:util';
import Papa from 'papaparse';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT, FORCE } from '../util/constants.js';
import { downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';

const exec = promisify(rawExec);
const execFile = promisify(callbackExecFile);


const DISCO_BATCH_URLS = [
  'https://zenodo.org/records/14159931/files/batch_1.tar.gz?download=1',
  'https://zenodo.org/records/14160154/files/batch_2.tar.gz?download=1',
  'https://zenodo.org/records/14160213/files/batch_3.tar.gz?download=1',
  'https://zenodo.org/records/14160221/files/batch_4.tar.gz?download=1',
  'https://zenodo.org/records/14160748/files/batch_5.tar.gz?download=1',
  'https://zenodo.org/records/14160802/files/batch_6.tar.gz?download=1',
  'https://zenodo.org/records/14166702/files/batch_7.tar.gz?download=1',
  'https://zenodo.org/records/15236185/files/batch_8.tar.gz?download=1',
  'https://zenodo.org/records/15236615/files/batch_9.tar.gz?download=1',
];

const DISCO_METADATA_URL = 'https://disco.bii.a-star.edu.sg/disco_v3_api/toolkit/getSampleMetadata';

const DISCO_BASE_URL = 'https://disco.bii.a-star.edu.sg/sample/';

const ORGAN_MAPPING = {
  // Add finalized mappings
  lung: 'UBERON:0002048',
  thymus: 'UBERON:0002370',
  eye: 'UBERON:0000970',
  uterus: 'UBERON:0000995',
  heart: 'UBERON:0000948',
  kidney: 'UBERON:0002113',
  trachea: 'UBERON:0003126',
  liver: 'UBERON:0002107',
  lymph_node: 'UBERON:0000029',
  spleen: 'UBERON:0002106',
  brain: 'UBERON:0000955',
  fallopian_tube: 'UBERON:0003889',
  ureter: 'UBERON:0000056',
  larynx: 'UBERON:0001737',
  pancreas: 'UBERON:0001264',
  spinal_cord: 'UBERON:0002240',
};

/** @implements {IDownloader} */
export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.cacheDir = getCacheDir(config);
    /** @type {string} */
    this.baseUrl = DISCO_BASE_URL;
    /** @type {string[]} */
    this.batchFiles = DISCO_BATCH_URLS.map((url) => {
      const batchName = url.match(/batch_\d+\.tar\.gz/)[0];
      return join(this.cacheDir, batchName);
    });
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.metadataFile = 'disco_metadata.tsv';
    /** @type {string} */
    this.metadataFilePath = join(this.cacheDir, this.metadataFile);
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(config, 'disco', this.extractScriptFile);
    /** @type {OrganMetadataCollection} */
    this.organMetadata = undefined;
  }

  async prepareDownload(datasets) {
    this.organMetadata = await OrganMetadataCollection.load(this.config);
    await downloadFile(this.metadataFilePath, DISCO_METADATA_URL, {
      overwrite: this.config.get(FORCE, false),
    });

    const metadataContent = fs.readFileSync(this.metadataFilePath, 'utf-8');
    const { data: records } = Papa.parse(metadataContent, {
      header: true,
      delimiter: '\t',
      skipEmptyLines: true,
    });

    // lookup table
    this.metadataLookup = {};
    for (const row of records) {
      this.metadataLookup['DISCO-' + row.sample_id] = row;
    }
  
    // Download all batch tar files
    await Promise.all(DISCO_BATCH_URLS.map((url) => this.downloadAndExtractBatch(url)));

    for (const dataset of datasets) {
      const sample_id = dataset.id.replace('DISCO-', '');
      dataset.dataset_id = `${this.baseUrl}${sample_id}`;
      dataset.consortium_name = 'DISCO';
      dataset.provider_name = 'DISCO';
      dataset.provider_uuid = 'bfe21f44-bf10-4371-8f4e-7f909f5a900c';
      dataset.dataset_link = `${this.baseUrl}${sample_id}`;
    }
  }

  async download(dataset) {
    // Parse the metadata

    // Find matching metadata row
    const matched = this.metadataLookup[dataset.id];
    if (!matched) {
      throw new Error(`No metadata found for dataset id: ${dataset.id}`);
    }

    // Assign cleaned metadata fields to dataset
    for (const [key, value] of Object.entries(matched)) {
      if (value !== undefined && value !== null && value.trim() !== '' && value.trim().toUpperCase() !== 'NA') {
        const targetKey = key === 'platform' ? 'dataset_technology' : key;
        dataset[targetKey] = value;
      }
    }

    // Resolve organ name from tissue using ORGAN_MAPPING
    const tissue = matched.tissue ?? '';
    const organCode = ORGAN_MAPPING[tissue] ?? '';
    dataset.organ = this.organMetadata.resolve(organCode);

    // Locate the .h5 file path for this sample
    const batchDirs = fs
      .readdirSync(this.cacheDir)
      .filter((name) => name.startsWith('batch_') && fs.statSync(join(this.cacheDir, name)).isDirectory());

    let h5FilePath = null;
    for (const dir of batchDirs) {
      const candidatePath = join(this.cacheDir, dir, `${matched.sample_id}.h5`);
      if (fs.existsSync(candidatePath)) {
        h5FilePath = candidatePath;
        break;
      }
    }

    if (!h5FilePath) {
      throw new Error(`Could not find .h5 file for ${dataset.id}`);
    }

    // Run Python extraction script
    const args = [
      this.extractScriptFilePath,
      '--metadata',
      this.metadataFilePath,
      '--dataset',
      h5FilePath,
      '--output',
      dataset.dataFilePath,
    ];
    const { stdout } = await execFile('python3', args);

    // Parse Python output and assign counts
    const counts = JSON.parse(stdout);
    dataset.dataset_cell_count = counts.cell_count;
    dataset.dataset_gene_count = counts.gene_count;

    // Enforce minimum cell count
    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cells. Cell count: ${dataset.dataset_cell_count}`);
    }
  }
  
  async downloadAndExtractBatch(url) {
    const batchName = url.match(/batch_\d+\.tar\.gz/)[0].replace('.tar.gz', '');
    const targetDir = join(this.cacheDir, batchName);

    if (existsSync(targetDir) && !this.config.get(FORCE, false)) {
      console.log(`Skipping already extracted batch: ${batchName}`);
      return;
    }

    fs.mkdirSync(targetDir, { recursive: true });

    const cmd = `curl -L "${url}" | tar -xz --strip-components=1 -C "${targetDir}"`;
    console.log(`Running command: ${cmd}`);
    await exec(cmd);
  }
}
