import { writeFile } from 'node:fs/promises';
import { getSummaryRef } from '../util/common.js';
import { concurrentMap } from '../util/concurrent-map.js';
import {
  ALGORITHMS,
  DEFAULT_MAX_CONCURRENCY,
  MAX_CONCURRENCY,
} from '../util/constants.js';
import { createSpec } from './spec.js';
import { getJobGeneratorRef } from './utils.js';

async function tryGenerateJobs(dataset) {
  const ref = getSummaryRef(dataset);
  try {
    const generator = getJobGeneratorRef(dataset);
    const metadata = await generator.createJob(dataset);
    ALGORITHMS.forEach((step) => {
      if (metadata[step] === false) {
        ref.setNotSupported(step);
      }
    });

    const spec = createSpec(metadata);
    const specString = JSON.stringify(spec, undefined, 2);
    await writeFile(dataset.jobFilePath, specString, { encoding: 'utf8' });
  } catch (error) {
    ALGORITHMS.forEach((step) => ref.setFailure(step, error.message ?? error));
  }
}

export async function generateJobs(datasets, config) {
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  await concurrentMap(datasets, tryGenerateJobs, { maxConcurrency });
}
