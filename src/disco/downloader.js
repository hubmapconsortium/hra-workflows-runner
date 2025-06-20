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

// const DEFAULT_DISCO_DOI = "https://disco.bii.a-star.edu.sg/";

const DISCO_BASE_URL = "https://disco.bii.a-star.edu.sg/"
// add more declrations as needed

const ORGAN_MAPPING = {

  // Add more mappings as needed
};

const execFile = promisify(callbackExecFile);

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
    for (let i = 0; i < DISCO_BATCH_URLS.length; i++) {
      const batchPath = this.batchFiles[i];

      // Check if file exists and we're not forcing
      if (existsSync(batchPath) && !this.config.get(FORCE, false)) {
        console.log(`Skipping existing file: ${batchPath}`);
        continue;
      }

      await downloadFile(
        batchPath, 
        DISCO_BATCH_URLS[i], 
        { overwrite: this.config.get(FORCE, false) }
      );
    }

    for (const dataset of datasets) {
      dataset.dataset_id = `${this.baseUrl}/sample/${dataset.id}`;
      dataset.publication = this.baseUrl;
      //dataset.publication_title = DISCO_PUBLICATION_NAME;
      //dataset.publication_lead_author = DISCO_PUBLICATION_LEAD_AUTHOR;
      dataset.consortium_name = 'DISCO';
      dataset.provider_name = 'DISCO';
      dataset.provider_uuid = 'disco-uuid'; // Update with actual UUID
      dataset.dataset_link = `${this.baseUrl}/sample/${dataset.id}`;
      dataset.dataset_technology = 'OTHER';
    }
  }

  async download(dataset) {
    const args = [
      this.extractScriptFilePath,
      '--metadata', this.metadataFilePath,
      '--dataset', dataset.id,
      '--output', dataset.dataFilePath,
      ...this.batchFiles
    ];

    console.log('Running Python with args:', ['python3', ...args].join(' '));

    const { stdout } = await execFile('python3', args);


    const metadata = JSON.parse(stdout);

    // Helper to assign non-empty fields
    const assignIfExists = (key) => {
      const value = metadata[key];
      if (value !== undefined && value !== null && value !== '') {
        dataset[key] = value;
      }
    };

    // Assign all fields from metadata
    assignIfExists("project_id");
    assignIfExists("sample_type");
    assignIfExists("tissue");

    // add logic for organ resolution (lowercasing and mapping)

    assignIfExists("anatomical_site");
    assignIfExists("disease");
    assignIfExists("platform");
    assignIfExists("age_group");
    assignIfExists("cell_sorting");
    assignIfExists("disease_subtype");
    assignIfExists("treatment");
    assignIfExists("time_point");
    assignIfExists("subject_id");
    assignIfExists("age");
    assignIfExists("gender");
    assignIfExists("race");
    assignIfExists("infection");
    assignIfExists("batch");
    assignIfExists("disease_stage");
    assignIfExists("genotype");
    assignIfExists("rna_source");
    assignIfExists("other_metadata");
    assignIfExists("disease_grade");
    assignIfExists("cell_number");
    assignIfExists("median_umi");
    assignIfExists("rds_md5");
    assignIfExists("rds_size");
    assignIfExists("source_cell_line");
    assignIfExists("source_tissue");
    assignIfExists("source_disease");
    assignIfExists("source_cell_type");
    assignIfExists("induced_cell_tissue");
    assignIfExists("collect_time");
    assignIfExists("sample_id");
    assignIfExists("dataset_cell_count");
    assignIfExists("dataset_gene_count");

    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cells. Cell count: ${dataset.dataset_cell_count}`);
    }
  }
}
