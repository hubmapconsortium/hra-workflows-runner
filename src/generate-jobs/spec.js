import { join } from 'node:path';

import { ALGORITHMS, DATA_FILE } from '../util/constants.js';
import { getCrosswalkingFilePath, getModelsDir } from '../util/paths.js';

const ALL_DISABLED_METADATA = ALGORITHMS.reduce((metadata, algorithm) => ({ ...metadata, [algorithm]: false }), {});

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
        tableLabelColumn: `${algorithm.charAt(0).toUpperCase() + algorithm.slice(1)}_Label`,
        tableClidColumn: 'CL_ID',
        tableMatchColumn: 'CL_Match',
      },
    }),
    summarize: {
      annotationMethod: algorithm,
      cellSource: metadata.cellSource,
    },
    directory: algorithm,
  };
}

export function createSpec(metadata, config, crosswalks) {
  const defaults = getAlgorithmDefaults(config);
  const algorithms = getEnabledAlgorithms(metadata);
  const algorithmSpecs = algorithms.map((algorithm) =>
    createAlgorithmSpec(config, algorithm, metadata, defaults, crosswalks[algorithm])
  );

  return {
    organ: metadata.organ,
    matrix: {
      class: 'File',
      path: DATA_FILE,
    },
    algorithms: algorithmSpecs,
  };
}

export function createSpecs(metadata, config, crosswalks) {
  const result = {};
  for (const algorithm of getEnabledAlgorithms(metadata)) {
    const newMetadata = {
      ...metadata,
      ...ALL_DISABLED_METADATA,
      [algorithm]: metadata[algorithm],
    };
    result[algorithm] = createSpec(newMetadata, config, crosswalks);
  }

  return result;
}
