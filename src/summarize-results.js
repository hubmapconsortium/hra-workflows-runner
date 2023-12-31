import { DatasetSummaries, DatasetSummary } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { Config } from './util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from './util/constants.js';
import { getAlgorithmReportFilePath, getDirForId, getSummariesFilePath } from './util/paths.js';

/**
 * Updates the summary status for each algorithm using the report files
 *
 * @param {DatasetSummary} item Summary item
 * @param {Config} config Configuration
 */
async function updateAlgorithmStatus(item, config) {
  const directory = getDirForId(item.id);
  item.errors = '';
  for (const algorithm of ALGORITHMS) {
    await readReport(item, algorithm, directory, config);
  }
}

/**
 * Reads a report file and updates the status for the algorithm
 *
 * @param {DatasetSummary} item Summary item
 * @param {string} algorithm Algorithm name
 * @param {string} directory Dataset directory name
 * @param {Config} config Configuration
 */
async function readReport(item, algorithm, directory, config) {
  const filePath = getAlgorithmReportFilePath(config, directory, algorithm);
  try {
    const { status, cause } = await loadJson(filePath);
    if (status === 'success') {
      item.setSuccess(algorithm);
    } else {
      item.setFailure(algorithm, cause);
    }
  } catch {
    // Ignore failures to load report
  }
}

async function main() {
  const config = getConfig();
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));

  await concurrentMap(Array.from(summaries.values()), (item) => updateAlgorithmStatus(item, config), {
    maxConcurrency,
  });

  await summaries.save(getSummariesFilePath(config));
}

main();
