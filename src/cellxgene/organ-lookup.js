import { readFileSync } from 'fs';
import Papa from 'papaparse';
import { checkFetchResponse } from '../util/fs.js';
import { logEvent } from '../util/logging.js';

/** crosswalk files parsed as jsonld */
const CROSSWALK_DATA = './crosswalking-tables/crosswalks.jsonld';
/** Endpoint to do tissue to organ query */
const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';

/**
 * Gets available organs
 *
 * @returns {Promise<string[]>}
 */
async function getOrgans() {
  const organs = new Set();
  const crosswalks = JSON.parse(readFileSync(CROSSWALK_DATA))['@graph'];
  for (const { organ_id: organ } of crosswalks) {
    if (organ?.trim().startsWith('UBERON')) {
      organs.add(organ);
    }
  }
  return Array.from(organs);
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
  console.log(query);
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
 * Creates a mapping from tissue id to organ id
 *
 * @param {string[]} ids Tissue ids
 * @returns {Promise<Map<string, string>>}
 */
export async function getOrganLookup(ids) {
  const organs = await logEvent('CellXGene:GetOrgans', () => getOrgans());
  const tissueIds = ids.filter((id) => !organs.includes(id));
  const organToOrganPairs = organs.map((organ) => [organ, organ]);
  const tissueToOrganPairs = await logEvent('CellXGene:GetOrganLookup', () => getLookup(organs, tissueIds));

  return new Map([...organToOrganPairs, ...tissueToOrganPairs]);
}
