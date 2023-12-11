import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { argv } from 'node:process';
import { DatasetSummaries } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { Config } from './util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY, SRC_DIR } from './util/constants.js';
import { getAlgorithmSummaryJsonLdFilePath, getDirForId, getOutputDir, getSummariesFilePath } from './util/paths.js';

/**
 * Loads a base context jsonld file
 *
 * @param {string} path Path to context file
 * @param {Config} config Configuration
 */
async function readContextFile(path, config) {
  const defaultPath = join(config.get(SRC_DIR), 'summary-context.jsonld');
  return await loadJson(path || defaultPath);
}

/**
 * Attempts to read the summary jsonld file for a specific dataset and algorithm
 *
 * @param {string} dir Dataset directory
 * @param {string} algorithm Algorithm name
 * @param {Config} config Configuration
 * @returns The summary file content or undefined if the file does not exist
 */
async function tryReadSummaryJsonLd(dir, algorithm, config) {
  const path = getAlgorithmSummaryJsonLdFilePath(config, dir, algorithm);
  try {
    return await loadJson(path);
  } catch {
    return undefined;
  }
}

/**
 * Reads all summary jsonld files for a dataset
 *
 * @param {DatasetSummaries} item Summary item
 * @param {Config} config Configuration
 */
async function readSummaryJsonLd(item, config) {
  const directory = getDirForId(item.id);
  return await concurrentMap(ALGORITHMS, (algorithm) => tryReadSummaryJsonLd(directory, algorithm, config));
}

async function main(contextFilePath) {
  const config = getConfig();
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const context = await readContextFile(contextFilePath, config);
  const jsonlds = await concurrentMap(Array.from(summaries.values()), (item) => readSummaryJsonLd(item, config), {
    maxConcurrency,
  });
  const entries = jsonlds
    .flat()
    .filter((json) => !!json)
    .flatMap((json) => json['@graph']);

  const outputFile = join(getOutputDir(config), 'bulk-cell-summaries.jsonld');
  const content = JSON.stringify(
    {
      ...context,
      '@graph': entries,
    },
    undefined,
    2
  );

  await writeFile(outputFile, content);
}

main(argv[2]);
