import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { argv } from 'node:process';
import { DatasetSummaries } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import {
  ALGORITHMS,
  DEFAULT_MAX_CONCURRENCY,
  MAX_CONCURRENCY,
  SRC_DIR,
} from './util/constants.js';
import {
  getAlgorithmSummaryJsonLdFilePath,
  getDirForId,
  getOutputDir,
  getSummariesFilePath,
} from './util/paths.js';

async function readContextFile(path, config) {
  const defaultPath = join(config.get(SRC_DIR), 'summary-context.jsonld');
  return await loadJson(path || defaultPath);
}

async function tryReadSummaryJsonLd(dir, algorithm, config) {
  const path = getAlgorithmSummaryJsonLdFilePath(config, dir, algorithm);
  try {
    return await loadJson(path);
  } catch {
    return undefined;
  }
}

async function readSummaryJsonLd(item, config) {
  const directory = getDirForId(item.id);
  return await concurrentMap(ALGORITHMS, (algorithm) =>
    tryReadSummaryJsonLd(directory, algorithm, config)
  );
}

async function main(contextFilePath) {
  const config = getConfig();
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const context = await readContextFile(contextFilePath, config);
  const jsonlds = await concurrentMap(
    Array.from(summaries.values()),
    (item) => readSummaryJsonLd(item, config),
    { maxConcurrency }
  );
  const entries = jsonlds
    .flat()
    .filter((json) => !!json)
    .flatMap((json) => json['@graph']);

  const outputFile = join(getOutputDir(config), 'summary.jsonld');
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
