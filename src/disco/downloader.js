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
import { getOrganLookup } from '../util/organ-lookup.js';

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

// Link to the mapping google sheet: https://docs.google.com/spreadsheets/d/1EkWBKOL-_YiR41MBv16w4KZzLxZ-0pgFx-FMJRJ5QiQ/edit?gid=470141504#gid=470141504
const TISSUE_MAPPING = {
  pituitary_gland: 'UBERON:0000007',
  lymph_node: 'UBERON:0000029',
  head: 'UBERON:0000033',
  tendon: 'UBERON:0000043',
  ureter: 'UBERON:0000056',
  intestine: 'UBERON:0000160',
  blood: 'UBERON:0000178',
  breast: 'UBERON:0000310',
  scalp: 'UBERON:0000403',
  testis: 'UBERON:0000473',
  stomach: 'UBERON:0000945',
  heart: 'UBERON:0000948',
  brain: 'UBERON:0000955',
  eye: 'UBERON:0000970',
  neck: 'UBERON:0000974',
  pleura: 'UBERON:0000977',
  leg: 'UBERON:0000978',
  ovary: 'UBERON:0000992',
  uterus: 'UBERON:0000995',
  seminal_vesicle: 'UBERON:0000998',
  nerve: 'UBERON:0001021',
  esophagus: 'UBERON:0001043',
  hypopharynx: 'UBERON:0001051',
  parathyroid_gland: 'UBERON:0001132',
  caecum: 'UBERON:0001153',
  colon: 'UBERON:0001155',
  peritoneal_cavity: 'UBERON:0001179',
  pancreas: 'UBERON:0001264',
  endometrium: 'UBERON:0001295',
  myometrium: 'UBERON:0001296',
  epididymis: 'UBERON:0001301',
  cerebrospinal_fluid: 'UBERON:0001359',
  arm: 'UBERON:0001460',
  ear: 'UBERON:0001690',
  nail: 'UBERON:0001705',
  nasal_cavity: 'UBERON:0001707',
  tongue: 'UBERON:0001723',
  nasopharynx: 'UBERON:0001728',
  oropharynx: 'UBERON:0001729',
  submandibular_gland: 'UBERON:0001736',
  larynx: 'UBERON:0001737',
  gingiva: 'UBERON:0001828',
  parotid_gland: 'UBERON:0001831',
  blood_vessel: 'UBERON:0001981',
  placenta: 'UBERON:0001987',
  thyroid_gland: 'UBERON:0002046',
  lung: 'UBERON:0002048',
  spleen: 'UBERON:0002106',
  liver: 'UBERON:0002107',
  gallbladder: 'UBERON:0002110',
  kidney: 'UBERON:0002113',
  duodenum: 'UBERON:0002114',
  jejunum: 'UBERON:0002115',
  ileum: 'UBERON:0002116',
  bronchiole: 'UBERON:0002186',
  spinal_cord: 'UBERON:0002240',
  brainstem: 'UBERON:0002298',
  umbilical_cord: 'UBERON:0002331',
  peritoneum: 'UBERON:0002358',
  adrenal_gland: 'UBERON:0002369',
  thymus: 'UBERON:0002370',
  bone_marrow: 'UBERON:0002371',
  tonsil: 'UBERON:0002372',
  bile_duct: 'UBERON:0002394',
  trachea: 'UBERON:0003126',
  omentum: 'UBERON:0003688',
  abdominal_wall: 'UBERON:0003697',
  fallopian_tube: 'UBERON:0003889',
  chest_wall: 'UBERON:0016435',
  umbilical_cord_blood: 'UBERON_0012168',
  gonad: 'UBERON:0000991',
  oral_cavity: 'UBERON:0000167',
  mucosa: 'UBERON:0000344',
  embryo: 'UBERON:0000922',
  urine: 'UBERON:0001088',
  mesonephros: 'UBERON:0000080',
  decidua: 'UBERON:0002450',
  nose: 'UBERON:0000004',
  dorsal_root_ganglion: 'UBERON:0000044',
  dental_pulp: 'UBERON:0001754',
  nucleus_pulposus: 'UBERON:0002242',
  corpus_cavernosum: 'UBERON:0006609',
  vestibulocochlear_nerve: 'UBERON:0001648',
  umbilical_vein: 'UBERON:0002066',
  annulus_fibrosus: 'UBERON:0006444',
  sympathetic_ganglion: 'UBERON:0001806',
  parietal_pleura: 'UBERON:0002400',
  vestibular_nerve: 'UBERON:0003723',
  pleural_fluid: 'UBERON:0001087',
};

// can have buckets for organs, tissues(this can be mapped to organ lookup)


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

    // Resolve organ name from tissue using TISSUE_MAPPING
    const tissue = matched.tissue ?? '';
    const tissueKey = tissue.replace(/\s+/g, '_'); // Replace all spaces with underscores
    const tissueUberonId = TISSUE_MAPPING[tissueKey] ?? '';
    //dataset.organ = this.organMetadata.resolve(organCode);

    // Changed Step : Tissue UBERON ID → organ UBERON ID (dynamic lookup)
    if (tissueUberonId) {
      const organLookup = await getOrganLookup([tissueUberonId], this.config, 'DISCO');
      dataset.organ = organLookup.get(tissueUberonId) ?? tissue;
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
