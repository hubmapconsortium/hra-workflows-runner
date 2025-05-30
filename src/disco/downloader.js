import { execFile as callbackExecFile } from 'node:child_process';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT } from '../util/constants.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';

const DISCO_BASE_URL = 'DISCO_BASE_URL';
const DEFAULT_DISCO_BASE_URL = 'https://path/to/disco/data'; // Update with actual base URL
const DISCO_DATA_MAPPING = 'DISCO_DATA_MAPPING';
const DEFAULT_DISCO_DATA_MAPPING = 'path/to/mapping.json'; // Update with actual mapping path
const DISCO_DOI_URL = 'https://doi.org/your-doi'; // Update with actual DOI
const DISCO_PUBLICATION_NAME = 'DISCO Publication Name'; // Update with actual publication name
const DISCO_PUBLICATION_LEAD_AUTHOR = 'Lead Author Name'; // Update with actual author name

const ORGAN_MAPPING = {
  // Add organ mappings similar to GTEx implementation
  // This should map DISCO tissue names to UBERON IDs
};

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.baseUrl = config.get(DISCO_BASE_URL, DEFAULT_DISCO_BASE_URL);
    /** @type {string} */
    this.mappingPath = config.get(DISCO_DATA_MAPPING, DEFAULT_DISCO_DATA_MAPPING);
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
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
      '--dataset', dataset.id,
      '--mapping-file', this.mappingPath,
      '--data-dir', getCacheDir(this.config),
      '--output', dataset.dataFilePath
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
