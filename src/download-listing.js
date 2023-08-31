import { DatasetSummaries, DatasetSummary } from './dataset/summary.js';
import { getConfig, loadListing } from './util/common.js';
import { DATASET_LIST_URL, FORCE } from './util/constants.js';
import { downloadFile, ensureDirsExist } from './util/fs.js';
import {
  getListingFilePath,
  getOutputDir,
  getSummariesFilePath,
} from './util/paths.js';

async function downloadListing(config) {
  const forceDownload = config.get(FORCE, false);
  const url = config.get(DATASET_LIST_URL);
  await downloadFile(getListingFilePath(config), url, {
    overwrite: forceDownload,
  });
}

function createSummaries(listing) {
  const items = listing.map((row) => new DatasetSummary(row.id));
  const summaries = new DatasetSummaries();
  summaries.push(...items);
  return summaries;
}

async function main() {
  const config = getConfig();
  await ensureDirsExist(getOutputDir(config));
  await downloadListing(config);

  const listing = await loadListing(config);
  const summaries = createSummaries(listing);
  await summaries.save(getSummariesFilePath(config));
}

main();
