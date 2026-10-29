import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT } from '../util/constants.js';
import { ensureDirsExist } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';
import {
  ANATOMOGRAM_BASE_URL,
  EXPERIMENTS,
  cacheExperiments,
  getExperimentFilePath,
  getExperimentFromId,
} from './utils.js';

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(config, 'anatomogram', this.extractScriptFile);
    /** @type {OrganMetadataCollection} */
    this.organMetadata = undefined;
  }

  async prepareDownload(datasets) {
    const config = this.config;
    this.organMetadata = await OrganMetadataCollection.load(config);

    await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));
    const experiments = await cacheExperiments(config);

    const supported = [];
    for (const dataset of datasets) {
      const code = getExperimentFromId(dataset.id);
      if (!experiments.includes(code)) {
        console.warn(`Anatomogram: Skipping ${dataset.id} from an unknown or disabled experiment`);
        continue;
      }

      const experiment = EXPERIMENTS[code];
      const baseIri = `${ANATOMOGRAM_BASE_URL}${code}#`;
      dataset.dataset_id = `${baseIri}${dataset.id}`;
      dataset.publication = experiment.publication;
      dataset.publication_title = experiment.publication_title;
      dataset.publication_lead_author = experiment.publication_lead_author;
      dataset.consortium_name = 'EBI Single Cell Expression Atlas';
      dataset.provider_name = 'EBI Single Cell Expression Atlas';
      dataset.provider_uuid = '148b56a3-a2f5-4b34-842a-c9d1e7666813';
      dataset.dataset_link = `${ANATOMOGRAM_BASE_URL}${code}`;
      dataset.dataset_technology = 'OTHER';
      dataset.dataset_rna_source = experiment.dataset_rna_source;
      supported.push(dataset);
    }

    return supported;
  }

  async download(dataset) {
    const code = getExperimentFromId(dataset.id);
    await ensureDirsExist(dataset.dirPath);

    const { stdout } = await execFile(
      'python3',
      [
        this.extractScriptFilePath,
        getExperimentFilePath(code, this.config),
        '--experiment',
        code,
        '--dataset',
        dataset.id,
        '--output',
        dataset.dataFilePath,
      ],
      { maxBuffer: 16 * 1024 * 1024 },
    );

    const metadata = JSON.parse(stdout);
    const baseIri = `${ANATOMOGRAM_BASE_URL}${code}#`;

    dataset.organ_source = metadata.organ ?? '';
    dataset.organ = this.organMetadata.resolve(EXPERIMENTS[code].organ);
    dataset.organ_id = dataset.organ ? `http://purl.obolibrary.org/obo/UBERON_${dataset.organ.split(':')[1]}` : '';
    dataset.donor_sex = metadata.sex ?? '';
    dataset.donor_age = metadata.age ?? '';
    dataset.donor_race = metadata.race ?? '';
    dataset.donor_disease = metadata.disease === 'normal' ? 'healthy' : (metadata.disease ?? '');
    dataset.donor_id = `${baseIri}${metadata.donor_id}`;
    dataset.dataset_cell_count = metadata.cell_count;
    dataset.dataset_gene_count = metadata.gene_count;
    dataset.block_id = `${dataset.dataset_id}_TissueBlock`;
    dataset.rui_location = '';

    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cell. Cell count: ${dataset.dataset_cell_count}`);
    }
  }
}
