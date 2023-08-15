import { mkdir } from 'node:fs/promises';
import { join } from 'node:path';
import { env } from 'node:process';
import papaparse from 'papaparse';
import { createDownloader } from './downloaders/index.js';
import { openFileWriteStream } from './downloaders/utils.js';
import { enrichDatasets } from './enrichers/index.js';

const DATASET_LIST_URL = env['DATASET_LIST_URL'];
const DATASET_ID_COLUMN = env['DATASET_ID_COLUMN'];
const DATASET_TYPE_HINT_COLUMN = env['DATASET_LINK_COLUMN'];

const OUTPUT_DIR = env['OUTPUT_DIR'] || '';
const CACHE_DIR = env['CACHE_DIR'] || 'tmp';
const ERROR_FILE = join(OUTPUT_DIR, 'errors.json');

const MAX_CONCURRENCY = Number.parseInt(env['MAX_CONCURRENCY'] || 10);

function checkEnv() {
  const required = [
    'DATASET_LIST_URL',
    'DATASET_ID_COLUMN',
    'DATASET_LINK_COLUMN',
    'OUTPUT_DIR',
  ];
  const missing = required.filter((key) => !env[key]);
  if (missing.length > 0) {
    const missingStr = missing.join(', ');
    const message = 'Missing required environment variables: ' + missingStr;
    throw new Error(message);
  }
}

async function fetchDatasetList(url) {
  const resp = await fetch(url);
  if (!resp.ok) {
    throw new Error(`Failed to download dataset list from '${url}'`);
  }

  return papaparse.parse(await resp.text(), {
    header: true,
    skipEmptyLines: 'greedy',
  });
}

function extractDatasets(data) {
  return data.map((row) => ({
    id: row[DATASET_ID_COLUMN],
    typeHint: row[DATASET_TYPE_HINT_COLUMN],
  }));
}

class Queue {
  constructor(items, maxLength = Infinity) {
    this.items = items;
    this.maxLength = maxLength;
    this.index = 0;
  }

  next() {
    const { items, maxLength, index } = this;
    return index < maxLength ? items[this.index++] : undefined;
  }
}

async function downloadExecutor(queue, errors) {
  let dataset;
  while ((dataset = queue.next()) !== undefined) {
    const downloader = createDownloader(dataset, OUTPUT_DIR, CACHE_DIR);
    if (!downloader) {
      errors.push({ dataset: dataset.id, cause: 'Not supported' });
      continue;
    }

    try {
      await downloader.download();
    } catch (error) {
      errors.push({ dataset: dataset.id, cause: error.message });
    }
  }
}

export async function main() {
  const errors = [];
  try {
    checkEnv();
    await mkdir(OUTPUT_DIR, { recursive: true });
    await mkdir(CACHE_DIR, { recursive: true });

    const { data } = await fetchDatasetList(DATASET_LIST_URL);
    const datasets = await enrichDatasets(extractDatasets(data));
    const queue = new Queue(datasets);
    const executors = Array(MAX_CONCURRENCY, errors)
      .fill(0)
      .map(() => downloadExecutor(queue, errors));

    await Promise.all(executors);
  } catch (error) {
    errors.push(error.message);
  }

  const stream = await openFileWriteStream(ERROR_FILE, { overwrite: true });
  stream?.end(JSON.stringify(errors, undefined, 2));
}

main();
