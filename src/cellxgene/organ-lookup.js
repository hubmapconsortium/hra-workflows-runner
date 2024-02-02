import Papa from 'papaparse';
import { checkFetchResponse } from '../util/fs.js';
import { logEvent } from '../util/logging.js';

/** Endpoint to get available organs */
const REF_ORGANS_ENDPOINT = 'https://grlc.io/api-git/hubmapconsortium/ccf-grlc/subdir/ccf//ref-organ-terms';
/** Endpoint to do tissue to organ query */
const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';

/** Extra organs not returned by getOrgans() */
const EXTRA_ORGANS = {
  blood: 'UBERON:0000178',
  bone_marrow: 'UBERON:0002371',
  adipose_tissue: 'UBERON:0001013',
};

/**
 * Gets available organs
 *
 * @returns {Promise<string[]>}
 */
async function getOrgans() {
  const resp = await fetch(REF_ORGANS_ENDPOINT, {
    headers: {
      Accept: 'application/json',
    },
  });
  checkFetchResponse(resp, 'CellxGene organ lookup failed');

  const items = await resp.json();
  return items
    .map((item) => item.representation_of)
    .filter((id) => id !== undefined)
    .map((id) => id.slice(id.lastIndexOf('/') + 1))
    .map((id) => id.replace('_', ':'))
    .filter((id) => id.startsWith('UBERON'));
}

/**
 * Get tissue to organ lookup pairs
 *
 * @param {string[]} organs Available organs
 * @param {string[]} ids Tissue ids
 * @returns {Promise<[string, string][]>}
 */
async function getLookup(organs, ids) {
  const serializeIds = (ids) => ids.map((id) => `(${id})`).join(' ');
  const query = `
    PREFIX part_of: <http://purl.obolibrary.org/obo/BFO_0000050>
    PREFIX UBERON: <http://purl.obolibrary.org/obo/UBERON_>

    SELECT DISTINCT ?as ?organ
    FROM <http://reasoner.renci.org/redundant>
    WHERE {
      VALUES (?organ_id) { ${serializeIds(organs)} }
      VALUES (?as_id) { ${serializeIds(ids)} }
      {
        ?as_id part_of: ?organ_id .
      }
      UNION
      {
        ?organ_id part_of: ?as_id .
      }

      BIND(REPLACE(STR(?as_id), 'http://purl.obolibrary.org/obo/UBERON_', 'UBERON:') as ?as)
      BIND(REPLACE(STR(?organ_id), 'http://purl.obolibrary.org/obo/UBERON_', 'UBERON:') as ?organ)
    }
  `;
  const body = new URLSearchParams({ query });
  const resp = await fetch(SPARQL_ENDPOINT, {
    method: 'POST',
    headers: {
      Accept: 'text/csv',
      'Content-Type': 'application/x-www-form-urlencoded',
      'Content-Length': body.toString().length.toString(),
    },
    body,
  });
  checkFetchResponse(resp, 'CellxGene organ lookup failed');

  const text = await resp.text();
  const {
    data: [_, ...tissueToOrganPairs],
  } = Papa.parse(text, { skipEmptyLines: 'greedy' });

  return tissueToOrganPairs;
}

/**
 * Creates a mapping from tissue id to organ id
 *
 * @param {string[]} ids Tissue ids
 * @returns {Promise<Map<string, string>>}
 */
export async function getOrganLookup(ids) {
  const organs = await logEvent('CellXGene:GetOrgans', () => getOrgans());
  const tissueIds = ids.filter((id) => !organs.includes(id));
  const organToOrganPairs = organs.map((organ) => [organ, organ]);
  const extraOrganPairs = Object.values(EXTRA_ORGANS).map((organ) => [organ, organ]);
  const tissueToOrganPairs = await logEvent('CellXGene:GetOrganLookup', () => getLookup(organs, tissueIds));

  return new Map([...organToOrganPairs, ...extraOrganPairs, ...tissueToOrganPairs]);
}
