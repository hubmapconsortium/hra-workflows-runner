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

// export class CollectionMetadata {
//   constructor(id, assets, donorTissuePairs, tissueIds) {
//     /** @type {string} */
//     this.id = id;
//     /** @type {string[]} */
//     this.assets = assets;
//     /** @type {[string, string][]} */
//     this.donorTissuePairs = donorTissuePairs;
//     /** @type {Map<string, string>} */
//     this.tissueIds = tissueIds;
//   }

//   findDataset(id) {
//     const isPrimary = ({ is_primary_data: value }) => /primary/i.test(value);
//     const isNonDiseased = ({ disease: [{ label }] }) => /normal/i.test(label);
//     const isHuman = ({ organism: [{ label }] }) => /homo sapiens/i.test(label);
//     return this.raw.datasets.find(
//       (dataset) => dataset.id === id // &&
//       // isPrimary(dataset) &&
//       // isNonDiseased(dataset) &&
//       // isHuman(dataset)
//     );
//   }

//   findH5adAsset(id) {
//     const isH5ad = ({ filetype }) => /^h5ad$/i.test(filetype);
//     const { dataset_assets: assets = [] } = this.findDataset(id) ?? {};
//     return assets.find(isH5ad);
//   }

//   findTissueId(id, tissue) {
//     const { tissue: tissues = [] } = this.findDataset(id) ?? {};
//     const entry = tissues.find(({ label }) => label === tissue);
//     return entry ? { id: entry.ontology_term_id, tissue } : undefined;
//   }

//   static async download(url) {
//     const resp = await fetch(url, { method: 'GET' });
//     checkFetchResponse(resp, 'CellXGene collection download failed');
//     return new CollectionMetadata(await resp.json());
//   }
// }
