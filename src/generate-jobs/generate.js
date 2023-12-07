import { writeFile } from 'node:fs/promises';
import { getSummaryRef } from '../util/common.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { createSpecs } from './spec.js';
import { getJobGeneratorRef } from './utils.js';
import { UnknownOrganError } from '../util/errors.js';
import { getCrosswalkingFilePath } from '../util/paths.js';
import { fileExists } from '../util/fs.js';

async function crosswalkExists(config) {
  const result = {};
  await concurrentMap(ALGORITHMS, async (algorithm) => {
    const path = getCrosswalkingFilePath(config, algorithm);
    result[algorithm] = await fileExists(path);
  });

  return result;
}

async function tryGenerateJobs(dataset, config) {
  const ref = getSummaryRef(dataset);
  try {
    const generator = getJobGeneratorRef(dataset);
    const metadata = await generator.createJob(dataset);
    ALGORITHMS.forEach((step) => {
      if (metadata[step] === false) {
        ref.setNotSupported(step);
      }
    });

    const crosswalks = await crosswalkExists(config);
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

export async function generateJobs(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  await concurrentMap(datasets, (dataset) => tryGenerateJobs(dataset, config), {
    maxConcurrency,
  });
}
