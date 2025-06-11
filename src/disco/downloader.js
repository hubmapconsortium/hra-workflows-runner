import { execFile as callbackExecFile } from 'node:child_process';
import { createHash } from 'node:crypto';
import { FORCE } from '../util/constants.js';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { OrganMetadataCollection } from '../organ/metadata.js';
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

const DISCO_DOI_URL = "https://doi.org/10.5281/zenodo.14159931";
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
    this.baseUrl = config.get(DISCO_BASE_URL, DEFAULT_DISCO_BASE_URL); 
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
      await downloadFile(
        this.batchFiles[i], 
        DISCO_BATCH_URLS[i], 
        { overwrite: this.config.get(FORCE, false) }
      );
    }

    for (const dataset of datasets) {
      dataset.dataset_id = `${DISCO_DOI_URL}#${dataset.id}`;
      dataset.publication = DISCO_DOI_URL;
      dataset.publication_title = DISCO_PUBLICATION_NAME;
      dataset.publication_lead_author = DISCO_PUBLICATION_LEAD_AUTHOR;
      dataset.consortium_name = 'DISCO';
      dataset.provider_name = 'DISCO';
      dataset.provider_uuid = 'your-provider-uuid'; // Update with actual UUID
      dataset.dataset_link = `${this.baseUrl}/datasets/${dataset.id}`;
      dataset.dataset_technology = 'OTHER';
    }
  }

  async download(dataset) {
    // Execute Python script to extract the dataset using the mapping file
    const { stdout } = await execFile('python3', [
      this.extractScriptFilePath,
      '--metadata', this.metadataFilePath,
      '--dataset', dataset.id,
      '--output', dataset.dataFilePath,
      // '--data-dir', getCacheDir(this.config),
      ...this.batchFiles
      
    ]);

    // Parse metadata from stdout
    const cell_count_match = /cell_count:\s*(\d+)\s*\n/i.exec(stdout);
    dataset.dataset_cell_count = parseInt(cell_count_match?.[1]);

    const gene_count_match = /gene_count:\s*(\d+)\s*\n/i.exec(stdout);
    dataset.dataset_gene_count = parseInt(gene_count_match?.[1]);

    const tissue_match = /tissue:(.+)\n/i.exec(stdout);
    const tissue = tissue_match?.[1].trim().toLowerCase() ?? '';
    dataset.organ = this.organMetadata.resolve(ORGAN_MAPPING[tissue] ?? '');
    dataset.organ_id = dataset.organ ? `http://purl.obolibrary.org/obo/UBERON_${dataset.organ.split(':')[1]}` : '';

    const donor_id_match = /donor_id:(.+)\n/i.exec(stdout);
    dataset.donor_id = `${DISCO_DOI_URL}#${donor_id_match?.[1].trim()}` ?? '';

    const sex_match = /sex:(.+)\n/i.exec(stdout);
    dataset.donor_sex = sex_match?.[1].trim() ?? '';

    const age_match = /age:(.+)\n/i.exec(stdout);
    dataset.donor_age_bin = age_match?.[1].trim() ?? '';

    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cells. Cell count: ${dataset.dataset_cell_count}`);
    }
  }
}
