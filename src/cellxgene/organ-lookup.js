import Papa from 'papaparse';
import { OrganMetadataCollection } from '../organ/metadata.js';
import { Config } from '../util/config.js';
import { checkFetchResponse } from '../util/fs.js';
import { logEvent } from '../util/logging.js';

/** Endpoint to do tissue to organ query */
const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';

/**
 * Get tissue to organ lookup pairs
 *
 * @param {string[]} organs Available organs
 * @param {string[]} ids Tissue ids
 * @returns {Promise<[string, string][]>}
 */
async function getTissueOrganPairs(organs, ids) {
  if (organs.length === 0 || ids.length === 0) {
    return [];
  }

  const serializeIds = (ids) => ids.map((id) => `(${id})`).join(' ');
  const query = `
    PREFIX part_of: <http://purl.obolibrary.org/obo/BFO_0000050>
    PREFIX UBERON: <http://purl.obolibrary.org/obo/UBERON_>
    PREFIX CL: <http://purl.obolibrary.org/obo/CL_>

    SELECT DISTINCT ?as ?organ
    FROM <http://reasoner.renci.org/redundant>
    WHERE {
      VALUES (?organ_id) { ${serializeIds(organs)} }
      VALUES (?as_id) { ${serializeIds(ids)} }

      ?as_id part_of: ?organ_id .

      BIND(REPLACE(REPLACE(STR(?as_id), STR(UBERON:), 'UBERON:'), STR(CL:), 'CL:') as ?as)
      BIND(REPLACE(STR(?organ_id), STR(UBERON:), 'UBERON:') as ?organ)
    }
  `;

  const resp = await fetch(SPARQL_ENDPOINT, {
    method: 'POST',
    headers: {
      Accept: 'text/csv',
      'Content-Type': 'application/sparql-query',
    },
    body: query,
  });
  checkFetchResponse(resp, 'CellxGene organ lookup failed');

  const text = await resp.text();
  const {
    data: [_, ...tissueToOrganPairs],
  } = Papa.parse(text, { skipEmptyLines: 'greedy' });

  return tissueToOrganPairs;
}

/**
 * Adds lookup from each organ to itself
 *
 * @param {OrganMetadataCollection} metadata Organ metadata
 * @param {Map<string, string>} lookup Map to add organs lookups
 */
function addOrgansToLookup(metadata, lookup) {
  for (const organ of metadata.organs) {
    lookup.set(organ, organ);
  }
}

/**
 * Adds a lookup from each tissue to it's organ
 * Selects the most specific organ if a tissue maps to multiple organs
 *
 * @param {[string, string][]} pairs Tissue to organ pairs
 * @param {OrganMetadataCollection} metadata Organ metadata
 * @param {Map<string, string>} lookup Map to add tissue to organ lookups
 */
function addTissueOrganPairsToLookup(pairs, metadata, lookup) {
  for (const [tissue, organ] of pairs) {
    const existing = lookup.get(tissue);
    const selectedOrgan = existing ? metadata.selectBySpecificity(organ, existing) : organ;
    lookup.set(tissue, selectedOrgan);
  }
}

/**
 * Resolves organs in the lookup mapping
 *
 * @param {OrganMetadataCollection} metadata Organ metadata
 * @param {Map<string, string>} lookup Map to add tissue to organ lookups
 */
function resolveOrgansInLookup(metadata, lookup) {
  for (const [key, organ] of lookup.entries()) {
    lookup.set(key, metadata.resolve(organ));
  }
}

/**
 * Creates a mapping from tissue id to organ id
 *
 * @param {string[]} ids Tissue ids
 * @param {Config} config Configuration
 * @returns {Promise<Map<string, string>>}
 */
export async function getOrganLookup(ids, config) {
  const metadata = await OrganMetadataCollection.load(config);
  const pairs = await logEvent('CellXGene:GetTissueOrganPairs', () =>
    getTissueOrganPairs(
      metadata.organs,
      ids.filter((id) => !metadata.has(id))
    )
  );
  const lookup = /** @type {Map<string, string>} */ new Map();
  addOrgansToLookup(metadata, lookup);
  addTissueOrganPairsToLookup(pairs, metadata, lookup);
  resolveOrgansInLookup(metadata, lookup);

  return lookup;
}
