import { ALGORITHMS, DATA_FILE } from '../util/constants.js';

const ALGORITHM_CELL_COLUMN = {
  azimuth: 'azimuth_label',
  celltypist: 'predicted_labels',
  popv: 'popv_prediction',
};

function createAlgorithmSpec(algorithm, metadata) {
  return {
    [algorithm]: metadata[algorithm] ?? {},
    extract: {
      annotationMethod: algorithm,
      cellLabelColumn: ALGORITHM_CELL_COLUMN[algorithm],
      geneLabelColumn: 'gene',
    },
  };
}

export function createSpec(metadata) {
  const algorithms = ALGORITHMS.filter(
    (algorithm) => metadata[algorithm] !== false
  ).map((algorithm) => createAlgorithmSpec(algorithm, metadata));

  return {
    organ: metadata.organ,
    matrix: {
      class: 'File',
      path: DATA_FILE,
    },
    preprocessing: {
      geneColumn: metadata.geneColumn,
    },
    algorithms,
  };
}
