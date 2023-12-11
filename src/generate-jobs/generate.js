import { writeFile } from 'node:fs/promises';
import { Dataset } from '../dataset/dataset.js';
import { getSummaryRef } from '../util/common.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { UnknownOrganError } from '../util/errors.js';
import { fileExists } from '../util/fs.js';
import { getCrosswalkingFilePath } from '../util/paths.js';
import { createSpecs } from './spec.js';
import { getJobGeneratorRef } from './utils.js';

/**
 * Checks for which algorithms a crosswalk file exists
 *
 * @param {Config} config Configuration
 * @returns {Promise<{ [algorithm: string]: boolean }>} Object mapping algorithm to whether it has a crosswalk file
 */
async function crosswalkExists(config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const result = {};

  await concurrentMap(
    ALGORITHMS,
    async (algorithm) => {
      const path = getCrosswalkingFilePath(config, algorithm);
      result[algorithm] = await fileExists(path);
    },
    { maxConcurrency }
  );

  return result;
}

/**
 * Generate jobs while handling errors and other failure conditions
 *
 * @param {Dataset[]} dataset Datasets to generate jobs for
 * @param {{ [algorithm: string]: boolean }} crosswalks Which crosswalks exist
 * @param {Config} config Configuration
 */
async function tryGenerateJobs(dataset, crosswalks, config) {
  const ref = getSummaryRef(dataset);
  try {
    const generator = getJobGeneratorRef(dataset);
    const metadata = await generator.createJob(dataset);
    ALGORITHMS.forEach((step) => {
      if (metadata[step] === false) {
        ref.setNotSupported(step);
      }
    });

    const specs = createSpecs(metadata, config, crosswalks);
    for (const algorithm in specs) {
      const spec = specs[algorithm];
      const specString = JSON.stringify(spec, undefined, 2);
      const filePath = dataset.jobFilePathWithSuffix(algorithm);
      await writeFile(filePath, specString, { encoding: 'utf8' });
    }
  } catch (error) {
    if (error instanceof UnknownOrganError) {
      ALGORITHMS.forEach((step) => ref.setNotSupported(step));
      ref.comments = error.message;
    } else {
      ALGORITHMS.forEach((step) => ref.setFailure(step, error.message ?? error));
    }
  }
}

/**
 * Generate job files for each dataset
 *
 * @param {Dataset[]} datasets Datasets to generate jobs for
 * @param {Config} config Configuration
 */
export async function generateJobs(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const crosswalks = await crosswalkExists(config);
  await concurrentMap(datasets, (dataset) => tryGenerateJobs(dataset, crosswalks, config), {
    maxConcurrency,
  });
}
