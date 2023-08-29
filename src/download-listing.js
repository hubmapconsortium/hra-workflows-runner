import { mkdir, readFile } from 'node:fs/promises';
import { join } from 'node:path';
import Papa from 'papaparse';
import { DatasetSummaries, DatasetSummary } from './dataset/summary.js';
import { Config } from './util/config.js';
import {
  DATASET_ID_COLUMN,
  DATASET_LINK_COLUMN,
  DATASET_LIST_URL,
  FORCE,
  LISTING_FILE,
  OUTPUT_DIR,
  SUMMARIES_FILE,
} from './util/constants.js';
import { defaultEnvReviver } from './util/env.js';
import { downloadFile } from './util/fs.js';

function createConfig() {
  return new Config()
    .loadEnv(defaultEnvReviver)
    .validate([
      OUTPUT_DIR,
      DATASET_LIST_URL,
      DATASET_ID_COLUMN,
      DATASET_LINK_COLUMN,
    ]);
}

async function downloadListing(config) {
  const forceDownload = config.get(FORCE, false);
  const listingPath = join(config.get(OUTPUT_DIR), LISTING_FILE);
  await downloadFile(listingPath, config.get(DATASET_LIST_URL), {
    overwrite: forceDownload,
  });
}

async function readListingIds(config) {
  const listingPath = join(config.get(OUTPUT_DIR), LISTING_FILE);
  const content = await readFile(listingPath, { encoding: 'utf8' });
  const listing = Papa.parse(content, { header: true }).data;
  const idColumn = config.get(DATASET_ID_COLUMN);
  return listing.map((row) => row[idColumn]);
}

async function createSummaries(ids, config) {
  const summariesPath = join(config.get(OUTPUT_DIR), SUMMARIES_FILE);
  const summaries = new DatasetSummaries();
  const items = ids.map((id) => new DatasetSummary(id));
  summaries.push(...items);
  await summaries.save(summariesPath);
}

async function main() {
  const config = createConfig();
  await mkdir(config.get(OUTPUT_DIR), { recursive: true });
  await downloadListing(config);

  const ids = await readListingIds(config);
  await createSummaries(ids, config);
}

main();
