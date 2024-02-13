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
  'https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad';
const GTEX_DOI_URL = 'https://doi.org/10.1126/science.abl4290';
const GTEX_PUBLICATION_NAME =
  'Single-nucleus cross-tissue molecular reference maps toward understanding disease gene function';
const GTEX_PUBLICATION_LEAD_AUTHOR = 'Gökcen Eraslan';
const GTEX_BLOCK_URL = 'https://gtexportal.org/home/tissue/';
const GTEX_PORTAL_LINK = 'https://gtexportal.org/home/singleCellOverviewPage';
const GTEX_EXTRACTION_SITES =
  'https://hubmapconsortium.github.io/hra-registrations/gtex-pan-eraslan-2022/rui_locations.jsonld';

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
  breast: 'UBERON:0001911',
  'esophagus mucosa': 'UBERON:0002469',
  'esophagus muscularis': 'UBERON:0004648',
  'skeletal muscle': 'UBERON:0001134',
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
    this.extractScriptFilePath = getSrcFilePath(config, 'gtex', this.extractScriptFile);
  }

  /**
   * Creates a lookup for extraction site
   *
   * @returns Lookup object
   */
  async fetchExtractionSiteLookup() {
    const req = await fetch(GTEX_EXTRACTION_SITES);
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
    await downloadFile(this.dataFilePath, this.dataUrl, {
      overwrite: this.config.get(FORCE, false),
    });

    for (const dataset of datasets) {
      dataset.dataset_id = `${GTEX_DOI_URL}#${dataset.id}`;
      dataset.publication = GTEX_DOI_URL;
      dataset.publication_title = GTEX_PUBLICATION_NAME;
      dataset.publication_lead_author = GTEX_PUBLICATION_LEAD_AUTHOR;
      dataset.consortium_name = 'GTEx';
      dataset.provider_name = 'GTEx';
      dataset.provider_uuid = '083882bb-6cc6-4c12-a205-eac37c1a2640';
      dataset.dataset_link = GTEX_PORTAL_LINK;
      dataset.dataset_technology = 'OTHER';
    }

    this.extractionSiteLookup = await this.fetchExtractionSiteLookup();
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
    dataset.organ = ORGAN_MAPPING[organ] ?? '';

    // Parse sex line. Format: `sex: X\n`
    const sex_match = /sex:(.+)\n/i.exec(stdout);
    dataset.donor_sex = sex_match?.[1].trim() ?? '';

    // Parse age line. Format: `age: X\n`
    const age_match = /age:(.+)\n/i.exec(stdout);
    dataset.donor_age_bin = age_match?.[1].trim() ?? '';

    // Parse donor_id line. Format: `donor_id: X\n`
    const donor_id_match = /donor_id:(.+)\n/i.exec(stdout);
    dataset.donor_id = `${GTEX_DOI_URL}#${donor_id_match?.[1].trim()}` ?? '';

    dataset.organ_id = dataset.organ ? `http://purl.obolibrary.org/obo/UBERON_${dataset.organ.split(':')[1]}` : '';

    // Parse tissue_site line. Format: `tissue_site: X\n`
    const tissue_site_match = /tissue_site:(.+)\n/i.exec(stdout);
    const tissueSite =
      `${GTEX_BLOCK_URL}${tissue_site_match?.[1]
        .trim()
        .replace(/[^a-zA-Z]+/g, '_')
        .replace(/_$/, '')}` ?? '';

    dataset.rui_location = this.extractionSiteLookup[tissueSite] ?? '';

    dataset.block_id = `${dataset.dataset_id}_TissueBlock`;
  }
}
