import { Dataset } from '../dataset/dataset.js';
import { Config } from '../util/config.js';
import { loadDatasetHandlers } from '../util/handler.js';
import { mapEntries } from '../util/iter.js';
import { getJobGeneratorRef, setJobGeneratorRef } from './utils.js';

/**
 * Attach job generators to each dataset
 *
 * @param {Dataset[]} datasets Datasets
 * @param {Config} config Configuration
 */
export async function attachGenerators(datasets, config) {
  const handlers = await loadDatasetHandlers(config);
  const generators = new Map(
    mapEntries(handlers, (handler) => new handler.JobGenerator(config))
  );

  for (const dataset of datasets) {
    const { handler } = dataset;
    const generator = generators.get(handler);

    setJobGeneratorRef(dataset, generator);
  }

  return datasets.filter(
    (dataset) => getJobGeneratorRef(dataset) !== undefined
  );
}
