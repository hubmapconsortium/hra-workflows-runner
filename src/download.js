import { writeFile } from 'fs/promises';
import Papa from 'papaparse';
import { DatasetSummaries } from './dataset/summary.js';
import { createDatasets } from './download/create-datasets.js';
import { download } from './download/download.js';
import { prepareDownloads } from './download/prepare-downloads.js';
import { getSummaryRef } from './download/utils.js';
import { getConfig, loadCsv } from './util/common.js';
import { DATASET_COLUMN_ID } from './util/constants.js';
import { ensureDirsExist } from './util/fs.js';
import { getCacheDir, getDataRepoDir, getListingFilePath, getSummariesFilePath } from './util/paths.js';

function filterSummaries(datasets) {
  const summaries = new DatasetSummaries();
  for (const dataset of datasets) {
    if (!dataset.scratch.exclude) {
      const summary = getSummaryRef(dataset);
      summary.organ = dataset.organ;
      summaries.push(summary);
    }
  }

  return summaries;
}

function filterListing(listing, summaries, config) {
  const idColumn = config.get(DATASET_COLUMN_ID);
  const summaryItems = Array.from(summaries.values());
  const ids = summaryItems.map((item) => item.id);
  const idSet = new Set(ids);
  const shouldIncludeRow = (item) => idSet.has(item[idColumn]);

  return listing.filter(shouldIncludeRow);
}

async function main() {
  const config = getConfig();
  await ensureDirsExist(getDataRepoDir(config), getCacheDir(config));

  const summariesFilePath = getSummariesFilePath(config);
  const summaries = await DatasetSummaries.load(summariesFilePath);

  const listingPath = getListingFilePath(config);
  const listing = await loadCsv(listingPath);

  const datasets = await createDatasets(summaries, config);
  const preparedDatasets = await prepareDownloads(datasets, config);
  await download(preparedDatasets, config);

  const newSummaries = filterSummaries(preparedDatasets);
  const newListing = filterListing(listing, newSummaries, config);
  const newListingCsv = Papa.unparse(newListing);

  await newSummaries.save(summariesFilePath);
  await writeFile(listingPath, newListingCsv);

  if (newListing.length < listing.length) {
    const diff = listing.length - newListing.length;
    console.info(`Excluded ${diff} datasets without any data`);
  }
}

main();
