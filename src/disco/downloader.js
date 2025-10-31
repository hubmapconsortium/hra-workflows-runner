import { execFile as callbackExecFile, exec as rawExec } from 'node:child_process';
import fs, { existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { promisify } from 'node:util';
import Papa from 'papaparse';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT, FORCE } from '../util/constants.js';
import { readCsv } from '../util/csv.js';
import { downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getOrganLookup } from '../util/organ-lookup.js';
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

const DISCO_METADATA_URL = 'https://www.immunesinglecell.com/disco_v3_api/toolkit/getSampleMetadata';

const DISCO_BASE_URL = 'https://www.immunesinglecell.com/sample/';

// Source Sheet can be downloaded from here: https://docs.google.com/spreadsheets/d/1EkWBKOL-_YiR41MBv16w4KZzLxZ-0pgFx-FMJRJ5QiQ/edit?gid=470141504#gid=470141504

async function loadTissueMappingFromCsv(csvPath) {
  if (!existsSync(csvPath)) {
    console.warn(`organ_mapping.csv not found at ${csvPath}; proceeding with empty mapping`);
    console.warn(`Please create organ_mapping.csv using the source sheet: https://docs.google.com/spreadsheets/d/1EkWBKOL-_YiR41MBv16w4KZzLxZ-0pgFx-FMJRJ5QiQ/edit?gid=470141504#gid=470141504`);
    return {};
  }
  /** @type {Record<string, string>} */
  const mapping = {};
  for await (const row of readCsv(csvPath)) {
    if (!row) continue;
    const rawLabel = (row.ontology_label ?? '').toString().trim();
    const rawId = (row.ontology_id ?? '').toString().trim();
    if (!rawLabel || !rawId) continue;
    const key = rawLabel
      .normalize('NFKC')
      .trim()
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, '_')
      .replace(/^_|_$/g, '');
    mapping[key] = rawId;
  }
  return mapping;
}

// Use CSV colocated with DISCO code (in same directory as this file)
const __dirname = dirname(fileURLToPath(import.meta.url));
const TISSUE_CSV = join(__dirname, 'organ_mapping.csv');
let TISSUE_MAPPING = null;

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
      dataset.donor_id = `${this.baseUrl}${sample_id}$Donor`;
      dataset.block_id = `${this.baseUrl}${sample_id}$TissueBlock`;
      dataset.rui_location = '';
      dataset.consortium_name = 'DISCO';
      dataset.provider_name = 'DISCO';
      dataset.provider_uuid = 'bfe21f44-bf10-4371-8f4e-7f909f5a900c';
      dataset.dataset_link = `${this.baseUrl}${sample_id}`;
      dataset.donor_link = `${this.baseUrl}${sample_id}`;
      dataset.block_link = `${this.baseUrl}${sample_id}`;
    }
  }

  async download(dataset) {
    // Parse the metadata & find matching metadata row
    const matched = this.metadataLookup[dataset.id];
    if (!matched) {
      throw new Error(`No metadata found for dataset id: ${dataset.id}`);
    }

    // Assign cleaned metadata fields to dataset
    for (let [key, value] of Object.entries(matched)) {
      if (value !== undefined && value !== null && value.trim() !== '' && value.trim().toUpperCase() !== 'NA') {
        let targetKey = key;
        switch (key) {
          case 'platform':
            targetKey = 'dataset_technology';
            break;
          case 'gender':
            targetKey = 'donor_sex';
            value = value === 'F' ? 'Female' : 'Male';
            break;
          case 'age':
            targetKey = 'donor_age';
            value = +value;
            break;
          case 'age_group':
            targetKey = 'donor_age_group';
            break;
          case 'race':
            targetKey = 'donor_race';
            break;
          case 'rna_source':
            targetKey = 'dataset_rna_source';
            break;
          case 'disease':
            targetKey = 'donor_disease';
            value = value === 'control' ? 'healthy' : value;
            break;
        }
        dataset[targetKey] = value;
      }
    }

    // Resolve organ name from tissue using TISSUE_MAPPING
    if (TISSUE_MAPPING === null) {
      TISSUE_MAPPING = await loadTissueMappingFromCsv(TISSUE_CSV);
    }
    const tissue = matched.tissue ?? '';
    const tissueKey = tissue
      .normalize('NFKC')
      .trim()
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, '_')
      .replace(/^_|_$/g, ''); // Normalize and replace all spaces with underscores
    const tissueUberonId = TISSUE_MAPPING[tissueKey] ?? '';
    if (tissueUberonId) {
      const organLookup = await getOrganLookup([tissueUberonId], this.config, 'DISCO');
      dataset.organ = organLookup.get(tissueUberonId) ?? tissueUberonId;
      dataset.organ_id = dataset.organ.replace('UBERON:', 'http://purl.obolibrary.org/obo/UBERON_');
    } else {
      throw new Error(`Could not determine organ for ${dataset.id}, tissue: ${tissue}`);
    }

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
