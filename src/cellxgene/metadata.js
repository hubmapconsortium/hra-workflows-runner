import { basename } from 'node:path';
import { checkFetchResponse } from '../util/fs.js';
import { groupBy } from '../util/iter.js';

/**
 * @typedef CollectionMetadata
 * @type {object}
 * @property {string} id
 * @property {{ id: string; dataset: string }[]} assets
 * @property {[string, string][]} donorTissuePairs
 * @property {Map<string, string>} tissueIdLookup
 */

/**
 * Parses metadata from an id
 *
 * @param {string} id Id to parse
 */
export function parseMetadataFromId(id) {
  const { pathname, hash } = new URL(id);
  const collection = basename(pathname);
  const decodedHash = decodeURIComponent(hash);
  const [donor, tissue] = decodedHash.slice(1).split('$');
  return { collection, donor, tissue };
}

export async function downloadCollectionMetadata(url) {
  const resp = await fetch(url, { method: 'GET' });
  checkFetchResponse(resp, 'CellXGene collection download failed');
  const raw = await resp.json();
  return parseCollectionMetadata(raw);
}

export function parseCollectionMetadata(raw) {
  const { id, datasets } = raw;
  const validDatasets = filterNonDiseasedHumanDatasets(datasets);
  const { primary, secondary } = partitionDatasetsByType(validDatasets);
  const selectedDatasets = selectDatasetUsingCellCount(primary, secondary);

  return /** @type {CollectionMetadata} */ ({
    id,
    assets: getAssets(selectedDatasets),
    donorTissuePairs: getDonorTissuePairs(selectedDatasets),
    tissueIdLookup: getTissueIdLookup(selectedDatasets),
  });
}

function filterNonDiseasedHumanDatasets(datasets) {
  const someLabelMatchesICase = (items, value) =>
    items.some(({ label }) => label.toLowerCase() === value);

  return datasets.filter(
    ({ disease, organism }) =>
      someLabelMatchesICase(disease, 'normal') &&
      someLabelMatchesICase(organism, 'homo sapiens')
  );
}

function partitionDatasetsByType(datasets) {
  const grouped = groupBy(datasets, ({ is_primary_data }) =>
    is_primary_data.toLowerCase()
  );

  return {
    primary: grouped.get('primary'),
    secondary: grouped.get('secondary'),
  };
}

function selectDatasetUsingCellCount(primary, secondary) {
  const sumCellCounts = (items) =>
    items.reduce((acc, { cell_count }) => acc + cell_count, 0);

  const primaryCount = sumCellCounts(primary);
  const secondaryCount = sumCellCounts(secondary);
  return primaryCount >= secondaryCount ? primary : secondary;
}

function getAssets(datasets) {
  return datasets
    .map(({ dataset_assets }) => dataset_assets)
    .flat()
    .filter(({ filetype }) => filetype.toLowerCase() === 'h5ad')
    .map(({ id, dataset_id }) => ({ id, dataset: dataset_id }));
}

function getDonorTissuePairs(datasets) {
  // TODO
  return [];
}

function getTissueIdLookup(datasets) {
  const items = datasets
    .map(({ tissue }) => tissue)
    .flat()
    .map(({ label, ontology_term_id }) => [
      label.toLowerCase(),
      ontology_term_id,
    ]);

  return new Map(items);
}
