import { execFile as callbackExecFile } from 'node:child_process';
import { createHash } from 'node:crypto';
import { FORCE } from '../util/constants.js';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { downloadFile } from '../util/fs.js';
import { existsSync } from 'node:fs';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT } from '../util/constants.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';
import fs from 'node:fs';
import Papa from 'papaparse' ;
import { spawn } from 'node:child_process';
import { pipeline } from 'node:stream';
import { createWriteStream } from 'node:fs';
import https from 'node:https';




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

const DISCO_BASE_URL = "https://disco.bii.a-star.edu.sg/sample/"


const ORGAN_MAPPING = {

  // Add more mappings as needed
};

const execFile = promisify(callbackExecFile);
const streamPipeline = promisify(pipeline);


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
    this.batchFiles = DISCO_BATCH_URLS.map(url => {
      const batchName = url.match(/batch_\d+\.tar\.gz/)[0];
      return join(this.cacheDir, batchName);
    });
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.metadataFile = "disco_metadata.tsv";
        /** @type {string} */
    this.metadataFilePath = join(this.cacheDir, this.metadataFile);
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(
      config,
      'disco',
      this.extractScriptFile
    );
    /** @type {OrganMetadataCollection} */
    this.organMetadata = undefined;
  }

  async prepareDownload(datasets) {

    this.organMetadata = await OrganMetadataCollection.load(this.config);
      await downloadFile(this.metadataFilePath, DISCO_METADATA_URL, {
        overwrite: this.config.get(FORCE, false)
      });

    // Download all batch tar files

    await Promise.all(DISCO_BATCH_URLS.map(async (url) =>  {

    // Extract batch name 
    const batchName = url.match(/batch_\d+\.tar\.gz/)[0].replace('.tar.gz', '');
    const expectedExtractedDir = join(this.cacheDir, batchName);

    // Skip if already extracted
    if (existsSync(expectedExtractedDir) && !this.config.get(FORCE, false)) {
      console.log(`Skipping already extracted batch: ${batchName}`);
      return;
    }

    console.log(`Streaming and extracting: ${url}`);

    await new Promise((resolve, reject) => {
      https.get(url, (res) => {
        if (res.statusCode !== 200) {
          return reject(new Error(`Request failed with status ${res.statusCode}`));
        }

        // tars are nested as batch_X/batch_X/ therefore stripping components
        const targetDir = join(this.cacheDir, batchName);
        fs.mkdirSync(targetDir, { recursive: true });  // ensure target dir exists
        
        const tar = spawn('tar', ['-xz', '--strip-components=1', '-C', targetDir]);
        res.pipe(tar.stdin);

        tar.on('close', (code) => {
          if (code === 0) {
            resolve();
          } else {
            reject(new Error(`tar exited with code ${code}`));
          }
        });

        tar.on('error', reject);
      }).on('error', reject);
    });
  }));


    for (const dataset of datasets) {
      dataset.dataset_id = `${this.baseUrl}${dataset.id}`;
      dataset.consortium_name = 'DISCO';
      dataset.provider_name = 'DISCO';
      dataset.provider_uuid = 'bfe21f44-bf10-4371-8f4e-7f909f5a900c'; // Update with actual UUID
      dataset.dataset_link = `${this.baseUrl}${dataset.id}`;
      dataset.dataset_technology = 'OTHER';
    }
  }

  async download(dataset) {
      // Step 1: Parse the metadata file
    const metadataContent = fs.readFileSync(this.metadataFilePath, 'utf-8');
    const { data: records, errors } = Papa.parse(metadataContent, {
      header: true,
      delimiter: '\t',
      skipEmptyLines: true
    });


      // Step 2: Find metadata row for current dataset
    const matched = records.find(row => row.sample_id === dataset.id);
    if (!matched) {
      throw new Error(`No metadata found for dataset id: ${dataset.id}`);
    }

  // Step 3: Assign only non-NA and non-empty fields to dataset
  for (const [key, value] of Object.entries(matched)) {
    if (
      value !== undefined &&
      value !== null &&
      value.trim() !== '' &&
      value.trim().toUpperCase() !== 'NA'
    ) {
      dataset[key] = value;
    }
  }

    // Step 3.5: Locate the .h5 file path for this sample
    const batchDirs = fs.readdirSync(this.cacheDir)
      .filter(name => name.startsWith('batch_') && fs.statSync(join(this.cacheDir, name)).isDirectory());

    let h5FilePath = null;
    for (const dir of batchDirs) {
      const candidatePath = join(this.cacheDir, dir, `${dataset.id}.h5`);
      if (fs.existsSync(candidatePath)) {
        h5FilePath = candidatePath;
        break;
      }
    }

    if (!h5FilePath) {
      throw new Error(`Could not find .h5 file for ${dataset.id}`);
    }


   // Step 4: Run Python extraction script
    const args = [
      this.extractScriptFilePath,
      '--metadata', this.metadataFilePath,
      '--dataset', h5FilePath,
      '--output', dataset.dataFilePath,
    ];

    console.log('Running Python with args:', ['python3', ...args].join(' '));

    const { stdout } = await execFile('python3', args);


    // Step 5: Parse Python output and assign counts
    const counts = JSON.parse(stdout);
    dataset.dataset_cell_count = counts.cell_count;
    dataset.dataset_gene_count = counts.gene_count;

    // Enforce minimum cell count
    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cells. Cell count: ${dataset.dataset_cell_count}`);
    }
  }
}
