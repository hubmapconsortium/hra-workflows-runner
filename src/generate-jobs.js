import { Status, Step } from './dataset/summary.js';
import { attachGenerators } from './generate-jobs/attach-generators.js';
import { generateJobs } from './generate-jobs/generate.js';
import { prepareJobs } from './generate-jobs/prepare-jobs.js';
import {
  getConfig,
  loadDatasets,
  loadSummaries,
  saveSummaries,
} from './util/common.js';

async function main() {
  const config = getConfig();
  const summaries = await loadSummaries(config);
  const downloadedItems = summaries.filterByStatus(
    Step.DOWNLOADED,
    Status.SUCCESS
  );

  const datasets = await loadDatasets(downloadedItems, config);
  const datasetsWithGenerator = await attachGenerators(datasets, config);
  const preparedDatasets = await prepareJobs(datasetsWithGenerator, config);
  await generateJobs(preparedDatasets, config);

  await saveSummaries(summaries, config);
}

main();
