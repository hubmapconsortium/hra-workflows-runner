import { checkFetchResponse } from '../util/fs.js';
import Papa from 'papaparse';

const SPARQL_ENDPOINT = 'https://ubergraph.apps.renci.org/sparql';
const ORGAN_IDS = [
  'UBERON:0002097',
  'UBERON:0000955',
  'UBERON:0000059',
  'UBERON:0002509',
  'UBERON:0004537',
  'UBERON:0004548',
  'UBERON:0004549',
  'UBERON:0001303',
  'UBERON:0000948',
  'UBERON:0001302',
  'UBERON:0004538',
  'UBERON:0004539',
  'UBERON:0002107',
  'UBERON:0001004',
  'UBERON:0001264',
  'UBERON:0001270',
  'UBERON:0000079',
  'UBERON:0002108',
  'UBERON:0002106',
  'UBERON:0002370',
  'UBERON:0001223',
  'UBERON:0001222',
  'UBERON:0001333',
  'UBERON:0001255',
  'UBERON:0002240',
  'UBERON:0001737',
  'UBERON:0000995',
  'UBERON:0001987',
  'UBERON:0002182',
  'UBERON:0003126',
  'UBERON:0000178',
  'UBERON:0002371',
  'UBERON:0000970',
  'UBERON:0002048',
  'UBERON:0000029',
  'UBERON:0001911',
  'UBERON:0002367',
];

function serializeIds(ids) {
  return ids.map((id) => `(${id})`).join(' ');
}

function createQuery(ids) {
  return `
    PREFIX part_of: <http://purl.obolibrary.org/obo/BFO_0000050>
    PREFIX UBERON: <http://purl.obolibrary.org/obo/UBERON_>

    SELECT DISTINCT ?as ?organ
    FROM <http://reasoner.renci.org/redundant>
    WHERE {
      VALUES (?organ_id) { ${serializeIds(ORGAN_IDS)} }
      VALUES (?as_id) { ${serializeIds(ids)} }
      ?as_id part_of: ?organ_id .

      BIND(REPLACE(STR(?as_id), 'http://purl.obolibrary.org/obo/UBERON_', 'UBERON:') as ?as)
      BIND(REPLACE(STR(?organ_id), 'http://purl.obolibrary.org/obo/UBERON_', 'UBERON:') as ?organ)
    }
  `;
}

export async function getOrganLookup(ids) {
  const body = new URLSearchParams({ query: createQuery(ids) });
  const headers = {
    Accept: 'text/csv',
    'Content-Type': 'application/x-www-form-urlencoded',
    'Content-Length': body.toString().length.toString(),
  };

  const resp = await fetch(SPARQL_ENDPOINT, {
    method: 'POST',
    headers,
    body,
  });
  checkFetchResponse(resp, 'CellxGene organ lookup failed');

  const text = await resp.text();
  const {
    data: [_, ...rows],
  } = Papa.parse(text, { skipEmptyLines: 'greedy' });

  return new Map(rows);
}
