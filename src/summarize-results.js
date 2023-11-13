import { DatasetSummaries, Status } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import {
  ALGORITHMS,
  DEFAULT_MAX_CONCURRENCY,
  MAX_CONCURRENCY,
} from './util/constants.js';
import {
  getAlgorithmReportFilePath,
  getDirForId,
  getSummariesFilePath,
  getDatasetFilePath,
} from './util/paths.js';

async function updateAlgorithmStatus(item, config) {
  const directory = getDirForId(item.id);
  item.errors = '';
  for (const algorithm of ALGORITHMS) {
    await readReport(item, algorithm, directory, config);
  }
}

async function updateOrgan(item, config) {
  const directory = getDirForId(item.id);
  await readOrganFromDatasetInfo(item, directory, config);
}

async function readOrganFromDatasetInfo(item, directory, config) {
  const filePath = getDatasetFilePath(config, directory);
  try {
    const { organ } = await loadJson(filePath);
    item.organ = organ;
  } catch {
    // Ignore failures to load dataset info
  }
}

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
  const downloadedItems = Array.from(summaries.values()).filter(
    (item) => item.downloaded === Status.SUCCESS
  );

  await concurrentMap(
    downloadedItems,
    (item) => updateOrgan(item, config),
    (item) => updateAlgorithmStatus(item, config),
    { maxConcurrency }
  );

  await summaries.save(getSummariesFilePath(config));
}

main();
