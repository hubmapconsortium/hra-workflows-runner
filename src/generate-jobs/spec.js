import { join } from 'node:path';

import { ALGORITHMS, DATA_FILE } from '../util/constants.js';
import { getModelsDir } from '../util/paths.js';

const ALGORITHM_CELL_COLUMN = {
  azimuth: 'azimuth_label',
  celltypist: 'predicted_labels',
  popv: 'popv_prediction',
};

const ALL_DISABLED_METADATA = ALGORITHMS.reduce(
  (metadata, algorithm) => ({ ...metadata, [algorithm]: false }),
  {}
);

function getAlgorithmDefaults(config) {
  return {
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
  };
}

function getEnabledAlgorithms(metadata) {
  return ALGORITHMS.filter((algorithm) => metadata[algorithm] !== false);
}

function createAlgorithmSpec(algorithm, metadata, defaults) {
  return {
    [algorithm]: {
      ...defaults[algorithm],
      ...metadata[algorithm],
    },
    extract: {
      annotationMethod: algorithm,
      cellLabelColumn: ALGORITHM_CELL_COLUMN[algorithm],
      geneLabelColumn: 'gene',
    },
  };
}

export function createSpec(metadata, config) {
  const defaults = getAlgorithmDefaults(config);
  const algorithms = getEnabledAlgorithms(metadata);
  const algorithmSpecs = algorithms.map((algorithm) =>
    createAlgorithmSpec(algorithm, metadata, defaults)
  );

  return {
    organ: metadata.organ,
    matrix: {
      class: 'File',
      path: DATA_FILE,
    },
    preprocessing: {
      geneColumn: metadata.geneColumn,
    },
    algorithms: algorithmSpecs,
  };
}

export function createSpecs(metadata, config) {
  const result = {};
  for (const algorithm of getEnabledAlgorithms(metadata)) {
    const newMetadata = {
      ...metadata,
      ...ALL_DISABLED_METADATA,
      [algorithm]: metadata[algorithm],
    };
    result[algorithm] = createSpec(newMetadata, config);
  }

  return result;
}
