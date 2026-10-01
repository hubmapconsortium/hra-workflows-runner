import { join } from 'node:path';

import { Config } from '../util/config.js';
import { ALGORITHMS, DATA_FILE, FILTERED_DATA_FILE, OPT_IN_ALGORITHMS } from '../util/constants.js';
import { getCrosswalkingFilePath, getModelsDir } from '../util/paths.js';

/** Metadata where all algorithms are disabled by default */
const ALL_DISABLED_METADATA = ALGORITHMS.reduce(
  (metadata, algorithm) => ({ ...metadata, [algorithm]: false }),
  /** @type {import('../util/handler.js').JobMetadata} */ ({})
);

/**
 * Get default values for each algorithm
 *
 * @param {Config} config Configuration
 * @returns Default values
 */
function getAlgorithmDefaults(config) {
  return {
    qc: {},
    azimuth: {
      referenceDataDir: {
        class: 'Directory',
        path: join(getModelsDir(config), 'azimuth'),
      },
    },
    celltypist: {},
    popv: {
      referenceDataDir: {
        class: 'Directory',
        path: join(getModelsDir(config), 'popv/reference-data'),
      },
      modelsDir: {
        class: 'Directory',
        path: join(getModelsDir(config), 'popv/models'),
      },
    },
    frmatch: {
      referenceDataDir: {
        class: 'Directory',
        path: join(getModelsDir(config), 'frmatch/reference-data'),
      },
      modelsDir: {
        class: 'Directory',
        path: join(getModelsDir(config), 'frmatch/models'),
      },
    },
    'pan-human-azimuth': {},
    author: {},
  };
}

/**
 * Checks whether an algorithm is enabled in the metadata.
 * Opt-in algorithms must be explicitly enabled by the job generator.
 *
 * @param {import('../util/handler.js').JobMetadata} metadata Metadata
 * @param {string} algorithm Algorithm name
 * @returns {boolean} Whether the algorithm is enabled
 */
export function isAlgorithmEnabled(metadata, algorithm) {
  if (OPT_IN_ALGORITHMS.includes(algorithm)) {
    return !!metadata[algorithm];
  }
  return metadata[algorithm] !== false;
}

/**
 * Get the enabled algorithms in the metadata
 *
 * @param {import('../util/handler.js').JobMetadata} metadata Metadata
 * @returns Names of enabled algorithms
 */
function getEnabledAlgorithms(metadata) {
  return ALGORITHMS.filter((algorithm) => isAlgorithmEnabled(metadata, algorithm));
}

/**
 * Creates the algorithm specification for a job
 *
 * @param {Config} config Configuration
 * @param {string} algorithm Algorithm name
 * @param {import('../util/handler.js').JobMetadata} metadata Metadata
 * @param {{ [algorithm: string]: object }} defaults Algorithm default values
 * @param {boolean} crosswalkExists Whether to add crosswalking to the job
 */
function createAlgorithmSpec(config, algorithm, metadata, defaults, crosswalkExists) {
  return {
    [algorithm]: {
      ...defaults[algorithm],
      ...metadata[algorithm],
    },
    ...(crosswalkExists && {
      crosswalk: {
        table: {
          class: 'File',
          path: getCrosswalkingFilePath(config, algorithm),
        },
        tableLabelColumn: 'Annotation_Label',
        tableClidColumn: 'CL_ID',
        tableMatchColumn: 'CL_Match',
      },
    }),
    geneExpression: {
      geneExprCount: Number(config.get('TOP_GENE_COUNT', 200)),
    },
    summarize: {
      annotationMethod: algorithm,
      cellSource: metadata.cellSource,
    },
    directory: algorithm,
  };
}

/**
 * Creates a single job specification with all enabled algorithms
 *
 * @param {import('../util/handler.js').JobMetadata} metadata Metadata
 * @param {Config} config Configuration
 * @param {{ [algorithm: string]: boolean }} crosswalks Whether crosswalk is enabled for each algorithm
 * @param {boolean} [useFilteredMatrix] Whether to use the QC filtered matrix as input
 */
export function createSpec(metadata, config, crosswalks, useFilteredMatrix = metadata.qc !== undefined) {
  const defaults = getAlgorithmDefaults(config);
  const algorithms = getEnabledAlgorithms(metadata);
  const algorithmSpecs = algorithms.map((algorithm) =>
    createAlgorithmSpec(config, algorithm, metadata, defaults, crosswalks[algorithm])
  );

  return {
    organ: metadata.organ,
    matrix: {
      class: 'File',
      path: useFilteredMatrix ? FILTERED_DATA_FILE : DATA_FILE,
    },
    algorithms: algorithmSpecs,
  };
}

/**
 * Creates a job specification for each enabled algorithm individually
 *
 * @param {import('../util/handler.js').JobMetadata} metadata Metadata
 * @param {Config} config Configuration
 * @param {{ [algorithm: string]: boolean }} crosswalks Whether crosswalk is enabled for each algorithm
 */
export function createSpecs(metadata, config, crosswalks) {
  const result = /** @type {{ [algorithm: string]: ReturnType<typeof createSpec> }} */ ({});
  const qcEnabled = isAlgorithmEnabled(metadata, 'qc');
  for (const algorithm of getEnabledAlgorithms(metadata)) {
    const newMetadata = {
      ...metadata,
      ...ALL_DISABLED_METADATA,
      [algorithm]: metadata[algorithm],
    };
    result[algorithm] = createSpec(newMetadata, config, crosswalks, algorithm !== 'qc' && qcEnabled);
  }

  return result;
}
