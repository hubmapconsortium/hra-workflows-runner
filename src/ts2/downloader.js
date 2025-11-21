import { execFile as callbackExecFile } from 'node:child_process';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT } from '../util/constants.js';
import { ensureDirsExist } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getDataRepoDir, getSrcFilePath } from '../util/paths.js';
import { cacheCollections } from './utils.js';

const TS2_FIGSHARE_ID = 'TS2_FIGSHARE_ID';
const DEFAULT_TS2_FIGSHARE_ID = '27921984';

const TS2_DATA_DOI = 'https://doi.org/10.6084/m9.figshare.27921984';
const TS2_PUBLICATION_DOI = 'https://doi.org/10.1101/2024.12.03.626516';
const TS2_PUBLICATION_NAME = 'Tabula Sapiens v2';
const TS2_PUBLICATION_LEAD_AUTHOR = 'Angela Pisco';
const TS2_PORTAL_LINK = TS2_PUBLICATION_DOI;
const TS2_EXTRACTION_SITES =
  'https://hubmapconsortium.github.io/hra-registrations/hca-pan-tabula_sapiens-2022/rui_locations.jsonld';

const ORGAN_MAPPING = {
  Bladder: 'UBERON:0001255',
  Blood: 'UBERON:0000178',
  Bone_Marrow: 'UBERON:0002371',
  Ear: '',
  Eye: 'UBERON:0000970',
  Fat: 'UBERON:0001013', // adipose
  Heart: 'UBERON:0000948',
  Kidney: 'UBERON:0002113',
  Large_Intestine: 'UBERON:0000059',
  Liver: 'UBERON:0002107',
  Lung: 'UBERON:0002048',
  Lymph_Node: 'UBERON:0000029', // or mesenteric lymph node (UBERON:0002509)?
  Mammary: 'UBERON:0001911',
  Muscle: '',
  Ovary: 'UBERON:0000992',
  Pancreas: 'UBERON:0001264',
  Prostate: 'UBERON:0002367',
  Salivary_Gland: '',
  Skin: 'UBERON:0002097',
  Small_Intestine: 'UBERON:0002108',
  Spleen: 'UBERON:0002106',
  Stomach: '',
  Testis: '',
  Thymus: 'UBERON:0002370',
  Tongue: 'UBERON:0001723',
  Trachea: 'UBERON:0003126',
  Uterus: 'UBERON:0000995',
  Vasculature: 'UBERON:0004537',
};

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.figshareId = config.get(TS2_FIGSHARE_ID, DEFAULT_TS2_FIGSHARE_ID);
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(config, 'ts2', this.extractScriptFile);
    /** @type {OrganMetadataCollection} */
    this.organMetadata = undefined;
  }

  /**
   * Creates a lookup for extraction site
   *
   * @returns Lookup object
   */
  async fetchExtractionSiteLookup() {
    const req = await fetch(TS2_EXTRACTION_SITES);
    const donors = await req.json();
    const lookup = /** @type {{ [tissueSite: string]: string }} */ ({});
    for (const donor of donors['@graph']) {
      for (const block of donor.samples) {
        const tissueSite = block.link;
        const extractionSite = JSON.stringify(block.rui_location);
        lookup[tissueSite] = extractionSite;
      }
    }
    return lookup;
  }

  async prepareDownload(datasets) {
    const config = this.config;
    this.organMetadata = await OrganMetadataCollection.load(config);

    await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));
    this.collections = await cacheCollections(this.figshareId, config);

    for (const dataset of datasets) {
      dataset.dataset_id = `${TS2_DATA_DOI}#${dataset.id}`;
      dataset.publication = TS2_PUBLICATION_DOI;
      dataset.publication_title = TS2_PUBLICATION_NAME;
      dataset.publication_lead_author = TS2_PUBLICATION_LEAD_AUTHOR;
      dataset.consortium_name = 'HCA';
      dataset.provider_name = 'Tabula Sapiens';
      dataset.provider_uuid = 'd65cd08e-7d9b-11ee-b962-0242ac120002';
      dataset.dataset_link = TS2_PORTAL_LINK;
      dataset.dataset_technology = 'OTHER';
      dataset.dataset_rna_source = "cell";
    }

    this.extractionSiteLookup = await this.fetchExtractionSiteLookup();
  }

  async download(dataset) {
    // Dataset ID Example: TS2-TSP14_Heart_NA_SS2_B134530_D101541_Live
    const organ_source = dataset.id
      // Clean up abbrevieated organs in the dataset id
      .replace('_LI_', '_Large_Intestine_')
      .replace('_SI_', '_Small_Intestine_')
      .replace('_BM_', '_Bone_Marrow_')
      .replace(/^TS2\-TSP[0-9]+\_/, '').toLowerCase().replace(/\_/g, '');
    const dataFile = this.collections.find((file) => {
      // Filename Example: Heart_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad
      const organ = file.name.split('_TSP')[0].toLowerCase().replace(/\_/g, '');
      return organ_source.startsWith(organ);
    })?.name;
    const dataFilePath = join(getCacheDir(this.config), 'ts2', dataFile);
    await ensureDirsExist(dataset.dirPath);
    
    const { stdout } = await execFile('python3', [
      this.extractScriptFilePath,
      dataFilePath,
      '--dataset',
      dataset.id,
      '--output',
      dataset.dataFilePath,
    ]);

    const metadata = JSON.parse(stdout);

    dataset.organ = this.organMetadata.resolve(ORGAN_MAPPING[metadata.organ] ?? '');
    dataset.organ_id = dataset.organ ? `http://purl.obolibrary.org/obo/UBERON_${dataset.organ.split(':')[1]}` : '';
    dataset.donor_sex = metadata.sex;
    dataset.donor_age = metadata.age;
    dataset.donor_race = metadata.race;
    dataset.donor_id = `${TS2_DATA_DOI}#${metadata.donor_id}` ?? '';
    dataset.dataset_cell_count = metadata.cell_count;
    dataset.dataset_gene_count = metadata.gene_count;
    dataset.block_id = `${dataset.dataset_id}_TissueBlock`;
    dataset.rui_location = this.extractionSiteLookup[metadata.tissue_site] ?? '';

    const minCount = this.config.get(DATASET_MIN_CELL_COUNT, DEFAULT_DATASET_MIN_CELL_COUNT);
    if (dataset.dataset_cell_count < minCount) {
      throw new Error(`Dataset has fewer than ${minCount} cell. Cell count: ${dataset.dataset_cell_count}`);
    }

    // TS2 datasets are generated from whole cells
    dataset.dataset_rna_source = metadata.dataset_rna_source ?? 'cell';
  }
}
