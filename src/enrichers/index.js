import enrichHubmapDatasets from './hubmap.js';

/**
 * @type {(function(object[]): Promise<object[] | undefined> | object[] | undefined)[]}
 */
const ENRICHERS = [enrichHubmapDatasets];

/**
 * Enrich datasets with additional data
 * @param {object[]} datasets All datasets
 * @returns The filtered and enriched datasets
 */
export async function enrichDatasets(datasets) {
  for (const fun of ENRICHERS) {
    datasets = (await fun(datasets)) ?? datasets;
  }

  return datasets;
}
