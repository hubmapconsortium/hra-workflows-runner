import { readFile, writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { argv } from 'node:process';
import Papa from 'papaparse';
import { DatasetSummaries, Status, Step } from './dataset/summary.js';
import { getConfig } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from './util/constants.js';
import {
  getAlgorithmSummaryJsonLdFilePath,
  getDirForId,
  getOutputDir,
  getSummariesFilePath,
} from './util/paths.js';
import { ensureDirsExist } from './util/fs.js';

async function readSummaryJsonLd(id, algorithm, config) {
  const path = getAlgorithmSummaryJsonLdFilePath(
    config,
    getDirForId(id),
    algorithm
  );
  const data = await readFile(path, { encoding: 'utf8' });
  return JSON.parse(data);
}

function countCellTypes(jsonld) {
  const items = jsonld['@graph'][0].summary;
  const countsByType = {};
  for (const { cell_label: type, count } of items) {
    countsByType[type] ??= 0;
    countsByType[type] += count;
  }

  return countsByType;
}

async function main(algorithm) {
  const config = getConfig();
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const downloadedItems = summaries.filterByStatus(
    Step.DOWNLOADED,
    Status.SUCCESS
  );
  const counts = await concurrentMap(
    downloadedItems,
    async (item) =>
      countCellTypes(await readSummaryJsonLd(item.id, algorithm, config)),
    { maxConcurrency }
  );
  const fields = Object.keys(
    counts.reduce((acc, cur) => ({ ...acc, ...cur }), {})
  ).sort();
  const zeroCounts = fields.reduce(
    (acc, field) => ({ ...acc, [field]: 0 }),
    {}
  );
  const countsWithIds = counts.map((counts, index) => ({
    ...zeroCounts,
    ...counts,
    id: downloadedItems[index].id,
  }));
  const resultCsv = Papa.unparse({
    fields: ['id', ...fields],
    data: countsWithIds,
  });
  const outputDir = join(getOutputDir(config), algorithm);
  const outputFilePath = join(outputDir, 'counts.csv');

  await ensureDirsExist(outputDir);
  await writeFile(outputFilePath, resultCsv, { encoding: 'utf8' });
}

main(argv[2]);
