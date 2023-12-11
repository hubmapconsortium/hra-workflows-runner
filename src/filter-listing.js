import { writeFile } from 'fs/promises';
import Papa from 'papaparse';
import { DatasetSummaries, DatasetSummary, Status, Step } from './dataset/summary.js';
import {
  getConfig,
  loadCsv,
  loadSummaries,
  saveSummaries,
} from './util/common.js';
import { DATASET_COLUMN_ID } from './util/constants.js';
import { getListingFilePath } from './util/paths.js';

/**
 * Tests whether the download step failed due to there not existing any data for the dataset
 *
 * @param {DatasetSummary} summary Summary data
 * @returns true if the download step errored due to not having any data, false otherwise
 */
function missingData(summary) {
  const error = summary.getError(Step.DOWNLOADED);
  return /^No data/i.test(error);
}

async function main() {
  const config = getConfig();

  const summaries = await loadSummaries(config);
  const failedItems = summaries.filterByStatus(Step.DOWNLOADED, Status.FAILURE);
  const missingDataItems = failedItems.filter(missingData);
  const missingDataItemsSet = new Set(missingDataItems);
  const items = Array.from(summaries.values());
  const newItems = items.filter((item) => !missingDataItemsSet.has(item));
  const newSummaries = new DatasetSummaries();
  newSummaries.push(...newItems);

  const listingFilePath = getListingFilePath(config);
  const listing = await loadCsv(listingFilePath);
  const idColumn = config.get(DATASET_COLUMN_ID);
  const missingDataIds = missingDataItems.map((item) => item.id);
  const missingDataIdsSet = new Set(missingDataIds);
  const newListingItems = listing.filter(
    (item) => !missingDataIdsSet.has(item[idColumn])
  );
  const newListing = Papa.unparse(newListingItems);

  await saveSummaries(newSummaries, config);
  await writeFile(listingFilePath, newListing);

  console.info(`Removed ${missingDataItems.length} datasets without any data`);
}

main();
