import { DatasetSummaries } from './dataset/summary.js';
import { createDatasets } from './download/create-datasets.js';
import { download } from './download/download.js';
import { prepareDownloads } from './download/prepare-downloads.js';
import { getConfig } from './util/common.js';
import { ensureDirsExist } from './util/fs.js';
import {
  getCacheDir,
  getDataRepoDir,
  getSummariesFilePath,
} from './util/paths.js';

async function main() {
  const config = getConfig();
  await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));

  const summariesFilePath = getSummariesFilePath(config);
  const summaries = await DatasetSummaries.load(summariesFilePath);
  const datasets = await createDatasets(summaries, config);
  const preparedDatasets = await prepareDownloads(datasets, config);
  await download(preparedDatasets, config);

  await summaries.save(summariesFilePath);
}

main();
