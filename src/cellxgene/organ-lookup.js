import { join } from 'node:path';
import Papa from 'papaparse';
import { loadJson } from '../util/common.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, FORCE, MAX_CONCURRENCY } from '../util/constants.js';
import { checkFetchResponse, downloadFile, ensureDirsExist } from '../util/fs.js';
import { logEvent } from '../util/logging.js';
import { getOutputDir } from '../util/paths.js';

/** Endpoint to do tissue to organ query */
const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';
/** Template for organ metadata file urls */
const ORGAN_METADATA_URL_TEMPLATE =
  'https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/containers/{{algorithm}}/context/organ-metadata.json';

/**
 * Downloads organ metadata for an algorithm
 *
 * @param {string} algorithm Algorithm organ metadata to download
 * @param {Config} config Configuration
 * @returns {Promise<Record<string, any>>} Organ metadata
 */
async function downloadAlgorithmOrganMetadata(algorithm, config) {
  const url = ORGAN_METADATA_URL_TEMPLATE.replace('{{algorithm}}', algorithm);
  const dir = join(getOutputDir(config), 'organ-metadata');
  const file = join(dir, `${algorithm}.json`);

  await ensureDirsExist(dir);
  await downloadFile(file, url, { overwrite: config.get(FORCE, false) });
  return await loadJson(file);
}

/**
 * Gets available organs
 *
 * @param {Config} config Configuration
 * @returns {Promise<string[]>}
 */
export async function getOrgans(config) {
  const metadata = await concurrentMap(ALGORITHMS, (algorithm) => downloadAlgorithmOrganMetadata(algorithm, config), {
    maxConcurrency: config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY),
  });

  return Array.from(new Set(metadata.flatMap(Object.keys)));
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
 * @param {Config} config Configuration
 * @returns {Promise<Map<string, string>>}
 */
export async function getOrganLookup(ids, config) {
  const organs = await logEvent('CellXGene:GetOrgans', () => getOrgans(config));
  const tissueIds = ids.filter((id) => !organs.includes(id));
  const organToOrganPairs = organs.map((organ) => [organ, organ]);
  const tissueToOrganPairs = await logEvent('CellXGene:GetOrganLookup', () => getLookup(organs, tissueIds));

  return new Map([...organToOrganPairs, ...tissueToOrganPairs]);
}
