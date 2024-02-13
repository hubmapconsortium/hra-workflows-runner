import { copyFile } from 'node:fs/promises';
import { DatasetSummaries, DatasetSummary } from './dataset/summary.js';
import { getConfig, loadListing } from './util/common.js';
import { Config } from './util/config.js';
import { DATASET_COLUMN_ID, DATASET_LIST_URL, FORCE } from './util/constants.js';
import { downloadFile, ensureDirsExist, fileExists } from './util/fs.js';
import { getDatasetListFilePath, getListingFilePath, getOutputDir, getSummariesFilePath } from './util/paths.js';

/**
 * Downloads listing from a remote source
 *
 * @param {Config} config Configuration
 */
async function downloadListing(config) {
  const forceDownload = config.get(FORCE, false);
  const url = config.get(DATASET_LIST_URL);
  await downloadFile(getListingFilePath(config), url, {
    overwrite: forceDownload,
  });
}

/**
 * Creates summary objects from listing rows
 *
 * @param {any[]} listing Listing rows
 */
function createSummaries(listing) {
  const items = listing.map((row) => new DatasetSummary(row.id));
  const summaries = new DatasetSummaries();
  summaries.push(...items);
  return summaries;
}

async function mergeExistingSummary(summaries, config) {
  const path = getSummariesFilePath(config);
  let existingSummaries;

  try {
    existingSummaries = await DatasetSummaries.load(path);
  } catch {
    // Failing to load existing summaries usually means they don't exist - Do nothing
    return;
  }

  for (const existingSummary of existingSummaries.values()) {
    const summary = summaries.get(existingSummary.id);
    if (summary) {
      Object.assign(summary, existingSummary);
    }
  }
}

async function main() {
  const config = getConfig().validate([DATASET_COLUMN_ID]);
  await ensureDirsExist(getOutputDir(config));

  const listFilePath = getDatasetListFilePath(config);
  const listingFilePath = getListingFilePath(config);
  if (await fileExists(listFilePath)) {
    await copyFile(listFilePath, listingFilePath);
  } else {
    config.validate([DATASET_LIST_URL]);
    await downloadListing(config);
  }

  const listing = await loadListing(config);
  const summaries = createSummaries(listing);
  await mergeExistingSummary(summaries, config);
  await summaries.save(getSummariesFilePath(config));
}

main();
