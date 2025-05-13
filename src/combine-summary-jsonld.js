import { createWriteStream } from 'node:fs';
import { join } from 'node:path';
import { createGzip } from 'node:zlib';
import { DatasetSummaries } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { Config } from './util/config.js';
import { ALGORITHMS } from './util/constants.js';
import { getAlgorithmSummaryJsonLdFilePath, getDirForId, getOutputDir, getSummariesFilePath } from './util/paths.js';

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

async function main() {
  const config = getConfig();
  const topGeneCount = Number(config.get('FINAL_TOP_GENE_COUNT', 10));
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const outputFile = join(getOutputDir(config), 'sc-transcriptomics-cell-summaries.jsonl.gz');
  let output = createWriteStream(outputFile, { autoClose: true });
  if (outputFile.endsWith('.gz')) {
    const gzip = createGzip();
    gzip.pipe(output);
    output = gzip;
  }
  for (const item of summaries.values()) {
    const jsonld = await readSummaryJsonLd(item, config);
    for (const entries of jsonld.flat().filter((s) => !!s)) {
      for (const entry of entries['@graph']) {
        for (const cellInfo of entry?.summary ?? []) {
          cellInfo['gene_expr'] = cellInfo['gene_expr'].slice(0, topGeneCount);
        }
        const content = JSON.stringify(entry);
        if (!output.write(content + '\n')) {
          // Drain buffer periodically
          await new Promise((resolve) => {
            output.once('drain', resolve);
          });
        }
      }
    }
  }
  output.end();
}

main();
