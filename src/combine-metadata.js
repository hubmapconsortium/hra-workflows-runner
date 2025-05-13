import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import Papa from 'papaparse';
import { DatasetSummaries, DatasetSummary } from './dataset/summary.js';
import { getConfig, loadJson } from './util/common.js';
import { concurrentMap } from './util/concurrent-map.js';
import { Config } from './util/config.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from './util/constants.js';
import { getDatasetFilePath, getDirForId, getOutputDir, getSummariesFilePath } from './util/paths.js';

/** Columns in the metadata csv */
const METADATA_FIELDS = [
  'id',
  'handler',
  'organ',
  'uuid',
  'assay_type',

  'dataset_id',
  'dataset_link',
  'dataset_technology',
  'dataset_info',
  'dataset_cell_count',
  'dataset_gene_count',

  'publication',
  'publication_title',
  'publication_lead_author',

  'consortium_name',
  'provider_name',
  'provider_uuid',

  'donor_id',
  'donor_age',
  'donor_development_stage',
  'donor_sex',
  'donor_bmi',
  'donor_race',

  'organ_id',
  'block_id',
  'section_id',

  'rui_location',
];

/**
 * Attempts to load metadata from file for a summary item
 *
 * @param {DatasetSummary} item Summary item
 * @param {Config} config Configuration
 * @returns Metadata or undefined if not found
 */
async function readMetadata(item, config) {
  const directory = getDirForId(item.id);
  const metadataFile = getDatasetFilePath(config, directory);

  try {
    return await loadJson(metadataFile);
  } catch {
    return undefined;
  }
}

function stringifyFields(obj, fields) {
  const shouldStringify = (value) => typeof value === 'object' && value !== null;
  const result = {};
  for (const field of fields) {
    const value = obj[field];
    result[field] = shouldStringify(value) ? JSON.stringify(value) : value;
  }

  return result;
}

async function main() {
  const config = getConfig();
  const maxConcurrency = config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY);
  const summaries = await DatasetSummaries.load(getSummariesFilePath(config));
  const metadata = await concurrentMap(Array.from(summaries.values()), (item) => readMetadata(item, config), {
    maxConcurrency,
  });
  const items = metadata.filter((item) => !!item).map((item) => stringifyFields(item, METADATA_FIELDS));

  const outputFile = join(getOutputDir(config), 'sc-transcriptomics-dataset-metadata.csv');
  const content = Papa.unparse({
    data: items,
    fields: METADATA_FIELDS,
  });

  await writeFile(outputFile, content);
}

main();
