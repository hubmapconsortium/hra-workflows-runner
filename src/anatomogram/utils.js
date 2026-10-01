import { join } from 'node:path';
import { concurrentMap } from '../util/concurrent-map.js';
import { FORCE } from '../util/constants.js';
import { downloadFile, ensureDirsExist } from '../util/fs.js';
import { getCacheDir } from '../util/paths.js';

const ANATOMOGRAM_EXPERIMENTS = 'ANATOMOGRAM_EXPERIMENTS';

export const ANATOMOGRAM_BASE_URL = 'https://www.ebi.ac.uk/gxa/sc/experiments/';
const ANATOMOGRAM_FTP_URL = 'https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/';

/** Metadata for each supported Single Cell Expression Atlas experiment */
export const EXPERIMENTS = {
  'E-CURD-119': {
    organ: 'UBERON:0002113', // kidney
    publication: 'https://doi.org/10.1038/s41467-021-22368-w',
    publication_title:
      'Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney',
    publication_lead_author: 'Yoshiharu Muto',
    dataset_rna_source: 'nucleus',
  },
  'E-MTAB-10553': {
    organ: 'UBERON:0002107', // liver
    publication: 'https://doi.org/10.1038/s41598-021-98806-y',
    publication_title:
      'Single-cell and bulk transcriptomics of the liver reveals potential targets of NASH with fibrosis',
    publication_lead_author: 'Zhong-Yi Wang',
    dataset_rna_source: 'cell',
  },
  'E-GEOD-130148': {
    organ: 'UBERON:0002048', // lung
    publication: 'https://doi.org/10.1038/s41591-019-0468-5',
    publication_title: 'A cellular census of human lungs identifies novel cell states in health and in asthma',
    publication_lead_author: 'Felipe A. Vieira Braga',
    dataset_rna_source: 'cell',
  },
  'E-MTAB-5061': {
    organ: 'UBERON:0001264', // pancreas
    publication: 'https://doi.org/10.1016/j.cmet.2016.08.020',
    publication_title: 'Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes',
    publication_lead_author: 'Åsa Segerstolpe',
    dataset_rna_source: 'cell',
  },
};

/**
 * Get the experiments to process
 *
 * @param {import('../util/config.js').Config} config Configuration
 * @returns {string[]} Experiment accession codes
 */
export function getExperiments(config) {
  const experiments = config.get(ANATOMOGRAM_EXPERIMENTS, Object.keys(EXPERIMENTS).join(','));
  return experiments
    .split(/[\s,;]+/)
    .filter((code) => code !== '')
    .filter((code) => {
      if (!(code in EXPERIMENTS)) {
        console.warn(`Anatomogram: Unknown experiment ${code}`);
        return false;
      }
      return true;
    });
}

/**
 * Parses the experiment code from a dataset id
 *
 * @param {string} id Dataset id, i.e. ANATOMOGRAM-E-CURD-119-Healthy1
 * @returns {string | undefined} Experiment code
 */
export function getExperimentFromId(id) {
  return /^ANATOMOGRAM-(E-[A-Z]+-\d+)-/i.exec(id)?.[1];
}

export function getExperimentFilePath(code, config) {
  return join(getCacheDir(config), 'anatomogram', `${code}.project.h5ad`);
}

/**
 * Downloads the h5ad file for each experiment into the cache
 *
 * @param {import('../util/config.js').Config} config Configuration
 * @returns {Promise<string[]>} The cached experiments
 */
export async function cacheExperiments(config) {
  const experiments = getExperiments(config);
  await ensureDirsExist(join(getCacheDir(config), 'anatomogram'));

  await concurrentMap(experiments, (code) =>
    downloadFile(getExperimentFilePath(code, config), `${ANATOMOGRAM_FTP_URL}${code}/${code}.project.h5ad`, {
      overwrite: config.get(FORCE, false),
    }),
  );

  return experiments;
}
