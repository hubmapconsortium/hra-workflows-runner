import { basename } from 'node:path';
import { checkFetchResponse } from '../util/fs.js';
import { groupBy } from '../util/iter.js';

/**
 * @typedef CollectionMetadata
 * @type {object}
 * @property {string} id
 * @property {{ id: string; dataset: string }[]} assets
 * @property {{ donor_id: string; tissue: string }[]} donorTissuePairs
 * @property {Map<string, string>} tissueIdLookup
 * @property {string} publication
 * @property {string} publication_name
 * @property {string} publication_lead_author
 * @property {string[]} consortium_name
 * @property {string} provider_name
 */

/** Portal endpoint */
const CELLXGENE_PORTAL_ENDPOINT = 'https://cellxgene.cziscience.com/collections/';

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

/**
 * Downloads collection metadata
 *
 * @param {string} url Collection metadata url
 * @returns Parsed metadata
 */
export async function downloadCollectionMetadata(url) {
  const resp = await fetch(url, { method: 'GET' });
  checkFetchResponse(resp, 'CellXGene collection download failed');
  const raw = await resp.json();
  return parseCollectionMetadata(raw);
}

/**
 * Parse a raw collection metadata object
 *
 * @param {any} raw Raw metadata object
 * @returns Parsed metadata
 */
export function parseCollectionMetadata(raw) {
  const { id, datasets, name } = raw;
  const validDatasets = filterNonDiseasedHumanDatasets(datasets);
  const { primary, secondary } = partitionDatasetsByType(validDatasets);
  const selectedDatasets = selectDatasetUsingCellCount(primary, secondary);
  return /** @type {CollectionMetadata} */ ({
    id,
    assets: getAssets(selectedDatasets),
    donorTissuePairs: getDonorTissuePairs(selectedDatasets),
    tissueIdLookup: getTissueIdLookup(selectedDatasets),
    dataset_link: `${CELLXGENE_PORTAL_ENDPOINT}${id}`,
    dataset_technology: 'OTHER',
    publication: getPublicationDOI(raw),
    publication_title: name,
    publication_lead_author: getPublicationLeadAuthor(raw),
    consortium_name: 'CxG', // alternate: consortia from raw
    provider_name: 'CxG', // alternate: curator_name from raw
    provider_uuid: 'f6841a8a-cef2-4421-b632-fc34ff5c27d8',
    assay_type: getAssayType(selectedDatasets),
  });
}

/**
 * Filter datasets that are from humans and without disease
 *
 * @param {any[]} datasets Raw dataset metadata
 * @returns Filtered datasets
 */
function filterNonDiseasedHumanDatasets(datasets) {
  const someLabelMatchesICase = (items, value) => items.some(({ label }) => label.toLowerCase() === value);

  return datasets.filter(
    ({ disease, organism }) =>
      someLabelMatchesICase(disease, 'normal') && someLabelMatchesICase(organism, 'homo sapiens')
  );
}

/**
 * Groups datasets by whether they are primary or secondary datasets
 *
 * @param {any[]} datasets Datasets to partition
 * @returns An object with the partitions
 */
function partitionDatasetsByType(datasets) {
  const grouped = groupBy(datasets, ({ is_primary_data }) => is_primary_data.toLowerCase());

  return {
    primary: grouped.get('primary') ?? [],
    secondary: grouped.get('secondary') ?? [],
  };
}

/**
 * Select datasets based on which group has the highest total cell count
 *
 * @param {any[]} primary Primary datasets
 * @param {any[]} secondary Secondary datasets
 * @returns Selected datasets
 */
function selectDatasetUsingCellCount(primary, secondary) {
  const sumCellCounts = (items) => items.reduce((acc, { cell_count }) => acc + cell_count, 0);

  const primaryCount = sumCellCounts(primary);
  const secondaryCount = sumCellCounts(secondary);
  return primaryCount >= secondaryCount ? primary : secondary;
}

/**
 * Extract h5ad datasets from dataset metadata
 *
 * @param {any[]} datasets Datasets
 * @returns {{id: string, dataset: string}[]}
 */
function getAssets(datasets) {
  return datasets
    .map(({ dataset_assets }) => dataset_assets)
    .flat()
    .filter(({ filetype }) => filetype.toLowerCase() === 'h5ad')
    .map(({ id, dataset_id }) => ({ id, dataset: dataset_id }));
}

/**
 * Computes unique donor-tissue pairs available in the datasets
 *
 * @param {any[]} datasets Datasets
 * @returns {{donor_id: string, tissue: string}[]} Unique donor-tissue pairs
 */
function getDonorTissuePairs(datasets) {
  const pairsArray = datasets.flatMap((dataset) =>
    dataset.donor_id.flatMap((donor_id) =>
      dataset.tissue.map((tissue) => ({
        donor_id: donor_id,
        tissue: tissue.label,
      }))
    )
  );

  const uniquePairs = Array.from(new Set(pairsArray.map(JSON.stringify)), JSON.parse);

  return uniquePairs;
}

/**
 * Creates a lookup from tissue label to ontoloty id
 *
 * @param {any[]} datasets Datasets
 * @returns {Map<string, string>}
 */
function getTissueIdLookup(datasets) {
  const items = datasets
    .map(({ tissue }) => tissue)
    .flat()
    .map(({ label, ontology_term_id }) => [label.toLowerCase(), ontology_term_id]);

  return new Map(items);
}

/**
 * Finds the publication doi url
 *
 * @param {any} raw Raw metadata
 * @returns {string}
 */
function getPublicationDOI(raw) {
  const { links } = raw;
  const doiLink = links.find((link) => link.link_type === 'DOI');
  return doiLink?.link_url ?? '';
}

/**
 * Finds the lead author of the publication
 *
 * @param {any} raw Raw metadata
 */
function getPublicationLeadAuthor(raw) {
  try {
    const {
      publisher_metadata: { authors: publication_authors },
    } = raw;
    return `${publication_authors[0].given} ${publication_authors[0].family}` ?? '';
  } catch {
    return '';
  }
}

/**
 * Get the unique assay types as a json encoded string
 *
 * @param {any[]} datasets Datasets
 */
function getAssayType(datasets) {
  const assayTypesArray = datasets.map(({ assay }) => assay).flat();
  return JSON.stringify(Array.from(new Set(assayTypesArray.map(JSON.stringify)), JSON.parse));
}
