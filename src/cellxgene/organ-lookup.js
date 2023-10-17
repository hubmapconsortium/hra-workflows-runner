import Papa from 'papaparse';
import { checkFetchResponse } from '../util/fs.js';

const REF_ORGANS_ENDPOINT =
  'https://grlc.io/api-git/hubmapconsortium/ccf-grlc/subdir/ccf//ref-organ-terms';
const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';

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

export async function getOrganLookup(ids) {
  const organs = await getOrgans();
  const tissueIds = ids.filter((id) => !organs.includes(id));
  const organToOrganPairs = organs.map((organ) => [organ, organ]);
  const tissueToOrganPairs = await getLookup(organs, tissueIds);

  return new Map([...organToOrganPairs, ...tissueToOrganPairs]);
}
